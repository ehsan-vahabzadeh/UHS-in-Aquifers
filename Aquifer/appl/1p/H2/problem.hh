// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup TwoPNCTests
 * \brief Definition of a problem for the 2pnc saline aquifer gas storage problem.
 */

#ifndef DUMUX_2PNC_TEST_PROBLEM_HH
#define DUMUX_2PNC_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <numeric>
#include <dumux/common/metadata.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/discretization/evalgradients.hh>
#include <dumux/io/json/json.hpp>

#include <stdexcept>

//! Exception thrown when a well pressure safety limit is violated
struct PressureLimitException : public std::runtime_error
{
    int exitCode;
    PressureLimitException(const std::string& msg, int code = 3)
        : std::runtime_error(msg), exitCode(code) {}
};

namespace Dumux {

/*!
 * \ingroup TwoPNCTests
 * \brief Two-phase n-component saline aquifer gas storage problem.
 *
 *  Migration note (1pnc to 2pnc):
 *
 *  Primary variables (TwoPNC, p_w-S_n formulation):
 *    [0] pressureIdx  - wetting-phase (water) pressure  p_w
 *    [1] switchIdx    - nonwetting-phase saturation S_n (when both phases
 *                       present) OR dissolved mole fraction (single-phase)
 *    [2..N-1]         - additional component mole fractions in wetting phase
 *
 *  Initial fully-saturated setup:
 *    S_w = 1 - eps,  S_n = eps  (eps = 1e-6 for numerical safety)
 *    Dissolved gas species initialised to zero (pure brine)
 *    Pressure initialised hydrostatically from Pressure_TOP / PressureGWC
 *
 *  Two-phase BC flux partitioning (production):
 *    Total component flux = Rate * sum_alpha [ x_alpha^kappa * lambda_alpha / sum lambda_alpha ]
 *    where alpha in {liquid, gas}, lambda_alpha = kr_alpha / mu_alpha (mobility),
 *    x_alpha^kappa = mole fraction of component kappa in phase alpha.
 *    This correctly weights both phases by their mobility.
 */
template <class TypeTag>
class OnePTwoCTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Extrusion = Extrusion_t<GridGeometry>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    enum
        {
            dimWorld = GridView::dimensionworld
        };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, 1>;
    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    
    // primary variable indices
        enum
        {
            pressureIdx = Indices::pressureIdx,
            switchIdx = Indices::switchIdx,
        };

        // phase indices
        enum
        {
            gasPhaseIdx = FluidSystem::gasPhaseIdx,
            liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
            H2OIdx = FluidSystem::H2OIdx,
            N2Idx = FluidSystem::N2Idx,
            CO2Idx = FluidSystem::CO2Idx,
            H2Idx = FluidSystem::H2Idx,
            CH4Idx = FluidSystem::CH4Idx,
        };
        // equation indices
        enum
        {
            numComponents = FluidSystem::numComponents,
            contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
            contiN2EqIdx = Indices::conti0EqIdx + N2Idx,
            contiCO2EqIdx = Indices::conti0EqIdx + CO2Idx,
            contiH2EqIdx = Indices::conti0EqIdx + H2Idx,
            contiCH4EqIdx = Indices::conti0EqIdx + CH4Idx,
        };
        // number of fluid phases
        enum { numPhases = FluidSystem::numPhases };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using Dimension_vector = Dune::FieldVector<Scalar, dimWorld>;
public:
    OnePTwoCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), useNitscheTypeBc_(getParam<bool>("Problem.UseNitscheTypeBc", false))
    {
        if (getParam<bool>("Problem.EnableGravity", true))
                gravity_ = {0, -9.81};
        //initialize fluid system
        FluidSystem::init();
        name_ = getParam<std::string>("Problem.Name");
        CELLS_VEC = getParam<Dimension_vector>("Grid.Cells");
        Upper_Right = getParam<Dimension_vector>("Grid.UpperRight");
        Lower_Left = getParam<Dimension_vector>("Grid.LowerLeft");
        X_max = Upper_Right[0];
        X_min = Lower_Left[0];
        Y_max = Upper_Right[1];
        Y_min = Lower_Left[1];
        Delta_x = (X_max - X_min);
        Delta_x = Delta_x / CELLS_VEC[0];
        Delta_y = (Y_max - Y_min);
        Delta_y = Delta_y / CELLS_VEC[1];
        // Position0 = getParam<Dimension_vector>("Grid.Positions0");
        // Position1 = getParam<Dimension_vector>("Grid.Positions1");
        // CELLS_VECX = getParam<Scalar>("Grid.Cells0");
        // CELLS_VECY = getParam<Scalar>("Grid.Cells1");
        // X_max = Position0[1];
        // X_min = Position0[0];
        // Y_max = Position1[1];
        // Y_min = Position1[0];
        // Delta_x = (X_max - X_min);
        // Delta_x = Delta_x / CELLS_VECX;
        // Delta_y = (Y_max - Y_min);
        // Delta_y = Delta_y / CELLS_VECY;
        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        // Pressure safety cutoff parameters
        enablePressureCutoff_ = getParam<bool>("Safety.EnablePressureCutoff", true);
        injPressureMultiplier_ = getParam<Scalar>("Safety.InjectionPressureMultiplier", 1.7);
        minProdPressure_ = getParam<Scalar>("Safety.MinProductionPressure", 1e5);
        pressureInitRef_ = getParam<Scalar>("Initialization.PressureGWC", 60e5);
    }

    const std::string &name() const
    {
        return name_;
    }
    /*!
     * \name Problem parameters
     */
    // \{

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        // values.setAllNeumann();
        if(globalPos[0] > X_max - eps_)
            values.setAllDirichlet();
        else
        values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);
        return values;
    }
     /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the tpfa scheme).
     */
     template<bool useBox = isBox, std::enable_if_t<!useBox, int> = 0>
     NumEqVector neumann(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        // no-flow everywhere except at the right boundary
        NumEqVector values(0.0);
        const auto xMax = this->gridGeometry().bBoxMax()[0];
        const auto& ipGlobal = scvf.ipGlobal();
        const auto &globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        Scalar Time = this->time();
        NumEqVector injectionComposition_(0.0);
        DimWorldMatrix molardensity(0.0), MW(0.0), Mobility(0.0);
        static const Scalar cylcesDev            = getParam<double>("BoundaryConditions.CyclesDev");
        static const Scalar injectionDurationDev = getParam<double>("BoundaryConditions.InjectionDurationDev")*86400;
        static const Scalar idleDurationDev     = getParam<double>("BoundaryConditions.IdleDurationDev")*86400;
        static const Scalar developmentDuration = cylcesDev*(injectionDurationDev+idleDurationDev);
        static const Scalar idleDurationOp     = getParam<double>("BoundaryConditions.IdleDurationOp")*86400;
        static const Scalar injectionDurationOp = getParam<double>("BoundaryConditions.InjectionDurationOp")*86400;
        static const Scalar extractionDurationOp = getParam<double>("BoundaryConditions.ExtractionDurationOp")*86400;
        injectionComposition_[H2Idx]            = getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
        injectionComposition_[CO2Idx]            = 1 - getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
        double inj_rate_dev                     = getParam<double>("BoundaryConditions.InjectionRateDev");
        const int cycleNumberOp                 = std::floor((Time - developmentDuration)/(injectionDurationOp+extractionDurationOp + idleDurationOp));
        const Scalar localTimeInCycle           = Time - developmentDuration - cycleNumberOp*(injectionDurationOp + extractionDurationOp + idleDurationOp);
        double inj_rate_Op = getParam<double>("BoundaryConditions.InjectionRateOp", 0.0);
        double Extractionrate_ = getParam<double>("BoundaryConditions.ProductionRate", 0.0);
        double inj_interval = getParam<double>("BoundaryConditions.Well_Height", 0.0);


        if ((globalPos[0] < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time < developmentDuration))
        // if ((globalPos[1] > Y_max - eps_)  && (Time < localTimeInCycle))
        {
            const int cycleNumber = std::floor(Time/(injectionDurationDev+idleDurationDev));
            const Scalar localTimeInCycle = Time - cycleNumber*(injectionDurationDev+idleDurationDev);
            if(localTimeInCycle <= injectionDurationDev){
                // values[CO2Idx] = (1-injectionComposition_[CO2Idx]) * inj_rate;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    // values[compIdx] = injectionComposition_[compIdx] * inj_rate_dev / FluidSystem::MixingFluidSystem::molarMass(compIdx);
                    values[compIdx] = injectionComposition_[compIdx] * inj_rate_dev ;
                }   
            }  
        }
        if ((globalPos[0] < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time >= developmentDuration))
        // if ((globalPos[1] > Y_max - eps_) && (Time >= developmentDuration))
        {
            if(localTimeInCycle < injectionDurationOp){
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    // values[compIdx] = injectionComposition_[compIdx] * inj_rate_Op / FluidSystem::MixingFluidSystem::molarMass(compIdx); 
                    values[compIdx] = injectionComposition_[compIdx] * inj_rate_Op;    
                }
            }
            if (localTimeInCycle > injectionDurationOp + idleDurationOp ){
                // Two-phase production: weight component flux by phase mobilities
                const Scalar mobLiq = volVars.mobility(liquidPhaseIdx);
                const Scalar mobGas = volVars.mobility(gasPhaseIdx);
                const Scalar totalMob = mobLiq + mobGas;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    const Scalar X_LIQ = volVars.moleFraction(liquidPhaseIdx, compIdx);
                    const Scalar X_GAS = volVars.moleFraction(gasPhaseIdx, compIdx);
                    if (totalMob > 0.0)
                        values[compIdx] = Extractionrate_ * (X_LIQ * mobLiq / totalMob + X_GAS * mobGas / totalMob);
                    else
                        values[compIdx] = 0.0;
                }
            }

        }
        
        if ((globalPos[0] > X_max - eps_))
        { 
            // set a fixed pressure on the right side of the domain
            auto init_val = initial_(globalPos);
            static const Scalar dirichletPressure = init_val[0];
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];

            auto d = ipGlobal - element.geometry().center();
            d /= d.two_norm2();

            const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d);

            // Two-phase outflow: sum over both phases
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                Scalar compFlux = 0.0;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    auto upwindTerm = useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx);
                    upwindTerm *= volVars.mobility(phaseIdx);
                    const auto phaseFlux = -1.0 * upwindTerm * tij * (dirichletPressure - volVars.pressure(phaseIdx));
                    compFlux += phaseFlux * volVars.moleFraction(phaseIdx, compIdx);
                }
                values[compIdx] = compFlux;
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the box scheme).
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<bool useBox = isBox, std::enable_if_t<useBox, int> = 0>
    NumEqVector neumann(const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const ElementVolumeVariables &elemVolVars,
                            const ElementFluxVariablesCache &elemFluxVarsCache,
                            const SubControlVolumeFace &scvf) const
    {   
        NumEqVector values(0.0);
        // return values;

        const auto &globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        Scalar Time = this->time();
        NumEqVector injectionComposition_(0.0);
        DimWorldMatrix molardensity(0.0), MW(0.0), Mobility(0.0);
        static const Scalar cylcesDev            = getParam<double>("BoundaryConditions.CyclesDev");
        static const Scalar injectionDurationDev = getParam<double>("BoundaryConditions.InjectionDurationDev")*86400;
        static const Scalar idleDurationDev     = getParam<double>("BoundaryConditions.IdleDurationDev")*86400;
        static const Scalar developmentDuration = cylcesDev*(injectionDurationDev+idleDurationDev);
        static const Scalar idleDurationOp     = getParam<double>("BoundaryConditions.IdleDurationOp")*86400;
        static const Scalar injectionDurationOp = getParam<double>("BoundaryConditions.InjectionDurationOp")*86400;
        static const Scalar extractionDurationOp = getParam<double>("BoundaryConditions.ExtractionDurationOp")*86400;
        // injectionComposition_[H2Idx]            = getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
        // injectionComposition_[CO2Idx]            = 1 - getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
        double inj_rate_dev                     = getParam<double>("BoundaryConditions.InjectionRateDev");
        const int cycleNumberOp                 = std::floor((Time - developmentDuration)/(injectionDurationOp+extractionDurationOp + idleDurationOp));
        const Scalar localTimeInCycle           = Time - developmentDuration - cycleNumberOp*(injectionDurationOp + extractionDurationOp + idleDurationOp);
        double inj_rate_Op = getParam<double>("BoundaryConditions.InjectionRateOp", 0.0);
        double Extractionrate_ = getParam<double>("BoundaryConditions.ProductionRate", 0.0);
        double inj_interval = getParam<double>("BoundaryConditions.Well_Height", 0.0);
        auto CushionGasType_ = getParam<std::string>("BoundaryConditions.CushionGasType");
        if (CushionGasType_ == "CH4")
            injectionComposition_[CH4Idx] = 1.0;
        else if (CushionGasType_ == "CO2")
            injectionComposition_[CO2Idx] = 1.0;
        else if (CushionGasType_ == "H2")
            injectionComposition_[H2Idx] = 1.0;
        else if (CushionGasType_ == "N2")
            injectionComposition_[N2Idx] = 1.0;
        else
        DUNE_THROW(Dune::InvalidStateException, "Invalid Cushion Gas Type " << CushionGasType_);

        if ((globalPos[0] < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time < developmentDuration))
    // if ((globalPos[1] > Y_max - eps_)  && (Time < localTimeInCycle))
        {
            const int cycleNumber = std::floor(Time/(injectionDurationDev+idleDurationDev));
            const Scalar localTimeInCycle = Time - cycleNumber*(injectionDurationDev+idleDurationDev);
            // std::cout << "injectionDurationDev: " << injectionDurationDev << std::endl;
            if(localTimeInCycle <= injectionDurationDev){
                // values[CO2Idx] = (1-injectionComposition_[CO2Idx]) * inj_rate;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    // values[compIdx] = injectionComposition_[compIdx] * inj_rate_dev / FluidSystem::MixingFluidSystem::molarMass(compIdx);
                    values[compIdx] = injectionComposition_[compIdx] * inj_rate_dev ;
                }   
            }  
        }
        if ((globalPos[0] < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time >= developmentDuration))
        // if ((globalPos[1] > Y_max - eps_) && (Time >= developmentDuration))
        {
            if(localTimeInCycle < injectionDurationOp){
                values[H2Idx] = inj_rate_Op;  
            }
            if (localTimeInCycle > injectionDurationOp + idleDurationOp ){
                // Two-phase production: weight component flux by phase mobilities
                const Scalar mobLiq = volVars.mobility(liquidPhaseIdx);
                const Scalar mobGas = volVars.mobility(gasPhaseIdx);
                const Scalar totalMob = mobLiq + mobGas;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    const Scalar X_LIQ = volVars.moleFraction(liquidPhaseIdx, compIdx);
                    const Scalar X_GAS = volVars.moleFraction(gasPhaseIdx, compIdx);
                    if (totalMob > 0.0)
                        values[compIdx] = Extractionrate_ * (X_LIQ * mobLiq / totalMob + X_GAS * mobGas / totalMob);
                    else
                        values[compIdx] = 0.0;
                }
            }

        }
        if ((globalPos[0] > X_max - eps_))
        {               
            // set a fixed pressure on the right side of the domain
            Scalar dirichletPressure(0.0);
            auto init_val = initial_(globalPos);
            const auto &volVars = elemVolVars[scvf.insideScvIdx()];
            const auto &fluxVarsCache = elemFluxVarsCache[scvf];
            dirichletPressure = init_val[0];

            // Two-phase outflow: sum contributions from both phases
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                Scalar compFlux = 0.0;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    // evaluate the pressure gradient for this phase
                    GlobalPosition gradP(0.0);
                    for (const auto &scv : scvs(fvGeometry))
                    {
                        const auto xIp = scv.dofPosition()[0];
                        auto tmp = fluxVarsCache.gradN(scv.localDofIndex());
                        tmp *= xIp > X_max - eps_ ? dirichletPressure
                                                    : elemVolVars[scv].pressure(phaseIdx);
                        gradP += tmp;
                    }

                    auto phaseFlux = vtmv(scvf.unitOuterNormal(), volVars.permeability(), gradP);
                    phaseFlux *= -1.0 * volVars.mobility(phaseIdx);
                    phaseFlux *= useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx);
                    compFlux += phaseFlux * volVars.moleFraction(phaseIdx, compIdx);
                }
                values[compIdx] = compFlux;
            }
        }

        return values;
    }

   
   void postTimeStep(const SolutionVector &curSol,
                          const SolutionVector &pastSol,
                          const GridVariables &gridVariables,
                          const Scalar timeStepSize,
                          std::shared_ptr<const GridGeometry> fvGridGeometry)
    {   
        // Scalar num_x = CELLS_VEC[0],num_grid = CELLS_VEC[0] * CELLS_VEC[1];
        Scalar t, dt, volume(0.0), numberOfCells(0.0);
        NumEqVector flow1(0.0);
        NumEqVector inventory(0.0), inventoryPastSol(0.0), materialBalanceError(0.0), inventory_error(0.0);
        Scalar averageReservoirPressure(0.0);
        // Scalar X_ID = 0, Z_ID = 0, Grid_ID = 0;   
        t = this->time();
        Scalar Time = t;
        dt = this->timeStepSize();
        NumEqVector values_inj(0.0), values_prod(0.0), outflow(0.0), values_inj_dt(0.0), values_prod_dt(0.0);
        
        for (const auto &element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*fvGridGeometry);
            fvGeometry.bind(element);

            // auto fvGeometry = localView(this->gridGeometry());
            // fvGeometry.bindElement(element);
            // auto globalPos = element.geometry().center();
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);
            auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
            prevElemVolVars.bindElement(element, fvGeometry, pastSol);

            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            const auto elemSolPast = elementSolution(element, pastSol, this->gridGeometry());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
            elemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (auto &&scv : scvs(fvGeometry))
            {
                const FluidState &fs = elemVolVars[scv].fluidState();
                const FluidState &fsOld = prevElemVolVars[scv].fluidState();
                // Scalar extrusion_  = extrusionFactor_sub(element, scv, elemSol);
                Scalar extrusion_ = Extrusion::volume(fvGeometry, scv);
                // Scalar extrusion_ = elemVolVars[scv].extrusionFactor();
                VolumeVariables volVars, volVarsPast;
                volVars.update(elemSol, *this, element, scv);
                volVarsPast.update(elemSolPast, *this, element, scv);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    
                inventory[compIdx] += volVars.porosity() * extrusion_  * ((fs.saturation(gasPhaseIdx) * fs.molarDensity(gasPhaseIdx) * fs.moleFraction(gasPhaseIdx, compIdx) + fs.saturation(liquidPhaseIdx) * fs.molarDensity(liquidPhaseIdx) * fs.moleFraction(liquidPhaseIdx, compIdx)));
                inventoryPastSol[compIdx] += volVarsPast.porosity() * extrusion_  * ((fsOld.saturation(gasPhaseIdx) * fsOld.molarDensity(gasPhaseIdx) * fsOld.moleFraction(gasPhaseIdx, compIdx) + fsOld.saturation(liquidPhaseIdx) * fsOld.molarDensity(liquidPhaseIdx) * fsOld.moleFraction(liquidPhaseIdx, compIdx)));
                }
                averageReservoirPressure += fs.pressure(gasPhaseIdx) * scv.volume();
                numberOfCells += 1;
                volume += scv.volume();
            }

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            NumEqVector injectionComposition_(0.0);
            // injectionComposition_[H2Idx] = 1.0;
            // injectionComposition_[CO2Idx] = 0.00;
            DimWorldMatrix molardensity(0.0), MW(0.0),Mobility(0.0);
            static const Scalar cylcesDev            = getParam<double>("BoundaryConditions.CyclesDev");
            static const Scalar injectionDurationDev = getParam<double>("BoundaryConditions.InjectionDurationDev")*86400;
            static const Scalar idleDurationDev     = getParam<double>("BoundaryConditions.IdleDurationDev")*86400;
            static const Scalar developmentDuration = cylcesDev*(injectionDurationDev+idleDurationDev);
            static const Scalar idleDurationOp     = getParam<double>("BoundaryConditions.IdleDurationOp")*86400;
            static const Scalar injectionDurationOp = getParam<double>("BoundaryConditions.InjectionDurationOp")*86400;
            static const Scalar extractionDurationOp = getParam<double>("BoundaryConditions.ExtractionDurationOp")*86400;
            // injectionComposition_[H2Idx]            = getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
            // injectionComposition_[CO2Idx]            = 1 - getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
            double inj_rate_dev                     = getParam<double>("BoundaryConditions.InjectionRateDev");
            const int cycleNumberOp                 = std::floor((Time - developmentDuration)/(injectionDurationOp+extractionDurationOp + idleDurationOp));
            const Scalar localTimeInCycle           = Time - developmentDuration - cycleNumberOp*(injectionDurationOp + extractionDurationOp + idleDurationOp);
            double inj_rate_Op = getParam<double>("BoundaryConditions.InjectionRateOp", 0.0);
            double Extractionrate_ = getParam<double>("BoundaryConditions.ProductionRate", 0.0);
            double inj_interval = getParam<double>("BoundaryConditions.Well_Height", 0.0);
            auto CushionGasType_ = getParam<std::string>("BoundaryConditions.CushionGasType");
            if (CushionGasType_ == "CH4")
                injectionComposition_[CH4Idx] = 1.0;
            else if (CushionGasType_ == "CO2")
                injectionComposition_[CO2Idx] = 1.0;
            else if (CushionGasType_ == "H2")
                injectionComposition_[H2Idx] = 1.0;
            else if (CushionGasType_ == "N2")
                injectionComposition_[N2Idx] = 1.0;
            else
            DUNE_THROW(Dune::InvalidStateException, "Invalid Cushion Gas Type " << CushionGasType_);
            for (auto&& scvf : scvfs(fvGeometry) )
            {   
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                VolumeVariables volVars, volVarsPast;
                volVars.update(elemSol, *this, element, scv);
                volVarsPast.update(elemSolPast, *this, element, scv);
                Scalar extrusion_  = Extrusion::area(fvGeometry, scvf);
                // if (scvf.boundary() && scvf.area() == Delta_y/2 && scvf.corner(0)[0] == X_min && scvf.corner(1)[0] == X_min)
                if (scvf.boundary()  && scvf.area() == Delta_y/2 && scvf.center()[0] == X_min)
                {
                    // const auto& globalPos = scvf.ipGlobal()[0];
                    // if ((globalPos[0] - Delta_x/2 < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time < localTimeInCycle))
                    if ((scvf.center()[1] - Y_min >= Y_max - inj_interval) && (Time < developmentDuration))
                    {
                        const int cycleNumber = std::floor(Time/(injectionDurationDev+idleDurationDev));
                        const Scalar localTimeInCycle = Time - cycleNumber*(injectionDurationDev+idleDurationDev);
                        if(localTimeInCycle <= injectionDurationDev){
                            // values[CO2Idx] = (1-injectionComposition_[CO2Idx]) * inj_rate;
                            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                            {
                                
                                // values_inj[compIdx] += injectionComposition_[compIdx] * Delta_y* inj_rate_dev / FluidSystem::MixingFluidSystem::molarMass(compIdx);
                                values_inj[compIdx] += injectionComposition_[compIdx] * extrusion_ * inj_rate_dev;
                            }
                        }
                    }
                    // if ((globalPos[0] - Delta_x/2 < X_min + eps_) && (globalPos[1] - Y_min >= Y_max - inj_interval) && (Time >= developmentDuration))
                    if ((scvf.center()[1] - Y_min >= Y_max - inj_interval) && (Time >= developmentDuration))
                    {
                        if(localTimeInCycle < injectionDurationOp){

                            values_inj[H2Idx] +=  extrusion_ * inj_rate_Op;  

                        }
                        if (localTimeInCycle > injectionDurationOp + idleDurationOp ){
                            // Two-phase production accounting
                            const Scalar mobLiq = volVars.mobility(liquidPhaseIdx);
                            const Scalar mobGas = volVars.mobility(gasPhaseIdx);
                            const Scalar totalMob = mobLiq + mobGas;
                            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                            {
                                const Scalar X_LIQ = volVars.moleFraction(liquidPhaseIdx, compIdx);
                                const Scalar X_GAS = volVars.moleFraction(gasPhaseIdx, compIdx);
                                if (totalMob > 0.0)
                                    values_prod[compIdx] += extrusion_ * Extractionrate_ * (X_LIQ * mobLiq / totalMob + X_GAS * mobGas / totalMob);
                            }
                        }
                    }

                
                }
            }
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        }
        const auto &comm = this->gridGeometry().gridView().comm();
            // Scalar sum_inj = std::accumulate(values_inj.begin(), values_inj.end(), 0);
            // Scalar sum_prod = std::accumulate(values_prod.begin(), values_prod.end(), 0);

            materialBalanceError = inventoryPastSol - inventory - dt * values_inj + dt * values_prod;

            // materialBalanceError = inventory + t * values_inj + t * values_prod;
            inventory_error = inventoryPastSol - inventory;
            Scalar error_sum = 0, inventory_error_sum = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
             {
                error_sum = error_sum + materialBalanceError[compIdx];
                inventory_error_sum = inventory_error_sum + inventory_error[compIdx];
                values_inj_dt[compIdx] = values_inj[compIdx] * dt;
                values_prod_dt[compIdx] = values_prod[compIdx] * dt;
            }
            // std::cout << "The material balance error is: " << materialBalanceError[2] <<  std::endl;
            // std::cout << "The material balance error is: " << materialBalanceError[2] <<  std::endl;

            // std::cout << "---------------------------------------------info end------------------------------------------------------------" << std::endl;

            Dumux::MetaData::Collector collector;
            if (Dumux::MetaData::jsonFileExists(name()))
                Dumux::MetaData::readJsonFile(collector, name());
            collector["time"].push_back(t);
            collector["averageReservoirPressure"].push_back(comm.sum(averageReservoirPressure) / comm.sum(volume));
            for (auto compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                collector["inventory"][FluidSystem::componentName(compIdx)].push_back(comm.sum(inventory[compIdx]));
                collector["materialBalanceError"][FluidSystem::componentName(compIdx)].push_back(comm.sum(materialBalanceError[compIdx]));
                collector["InjectionValues"][FluidSystem::componentName(compIdx)].push_back(comm.sum(values_inj[compIdx]));
                collector["ProductionValues"][FluidSystem::componentName(compIdx)].push_back(comm.sum(values_prod[compIdx]));
                collector["InjectionValues_dt"][FluidSystem::componentName(compIdx)].push_back(comm.sum(values_inj_dt[compIdx]));
                collector["ProductionValues_dt"][FluidSystem::componentName(compIdx)].push_back(comm.sum(values_prod_dt[compIdx]));
            }
            
            Dumux::MetaData::writeJsonFile(collector, name());

            // ---- Pressure safety checks ----
            if (enablePressureCutoff_)
            {
                double inj_interval = getParam<double>("BoundaryConditions.Well_Height", 0.0);

                Scalar pInjSum = 0.0, pProdSum = 0.0;
                int nInj = 0, nProd = 0;

                for (const auto &elem : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    auto fvGeo = localView(*fvGridGeometry);
                    fvGeo.bindElement(elem);
                    auto evv = localView(gridVariables.curGridVolVars());
                    evv.bindElement(elem, fvGeo, curSol);

                    for (auto&& scv : scvs(fvGeo))
                    {
                        const auto& pos = scv.dofPosition();
                        if (pos[0] < X_min + Delta_x && pos[1] - Y_min >= Y_max - inj_interval)
                        {
                            const auto& vv = evv[scv];
                            Scalar p = vv.pressure(liquidPhaseIdx);
                            pInjSum += p;
                            pProdSum += p;
                            nInj++;
                            nProd++;
                        }
                    }
                }

                const auto &comm = this->gridGeometry().gridView().comm();
                const Scalar pInjSumGlobal = comm.sum(pInjSum);
                const Scalar pProdSumGlobal = comm.sum(pProdSum);
                const int nInjGlobal = comm.sum(nInj);
                const int nProdGlobal = comm.sum(nProd);

                Scalar pInjMonitored = (nInjGlobal > 0) ? pInjSumGlobal / nInjGlobal : 0.0;
                Scalar pProdMonitored = (nProdGlobal > 0) ? pProdSumGlobal / nProdGlobal : 0.0;

                const Scalar maxInjPressure = injPressureMultiplier_ * pressureInitRef_;

                auto writeFailure = [&](const std::string& reason)
                {
                    if (comm.rank() == 0)
                    {
                        Dumux::MetaData::Collector fc;
                        if (Dumux::MetaData::jsonFileExists(name()))
                            Dumux::MetaData::readJsonFile(fc, name());
                        fc["runStatus"] = "failed_pressure_limit";
                        fc["failureReason"] = reason;
                        fc["failureTime"] = t;
                        fc["pInjMonitored"] = pInjMonitored;
                        fc["pProdMonitored"] = pProdMonitored;
                        fc["pressureInitRef"] = pressureInitRef_;
                        Dumux::MetaData::writeJsonFile(fc, name());
                    }
                };

                if (pInjMonitored > maxInjPressure)
                {
                    std::cout << "\n*** PRESSURE SAFETY: Injection pressure "
                              << pInjMonitored << " Pa exceeds limit "
                              << maxInjPressure << " Pa (= "
                              << injPressureMultiplier_ << " * " << pressureInitRef_
                              << ") at time " << t << " s ***\n" << std::endl;
                    writeFailure("inj_pressure_exceeded");
                    throw PressureLimitException("Injection pressure exceeded safety limit", 3);
                }

                if (pProdMonitored < minProdPressure_ && pProdMonitored > 0.0)
                {
                    std::cout << "\n*** PRESSURE SAFETY: Production pressure "
                              << pProdMonitored << " Pa below minimum "
                              << minProdPressure_ << " Pa at time " << t << " s ***\n" << std::endl;
                    writeFailure("prod_pressure_too_low");
                    throw PressureLimitException("Production pressure below safety limit", 4);
                }
            }
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilated per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions (mole/(m^3*s) or kg/(m^3*s)).
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}
    Scalar time() const
        {
            return timeLoop_->time();
        }

        void setTimeLoop(TimeLoopPtr timeLoop)
        {
            timeLoop_ = timeLoop;
            // timeLoop_->setCheckPoint({86400, 3 * 2.592e+6});
        }

        Scalar timeStepSize() const
        {
            return timeLoop_->timeStepSize();
        }
        const GlobalPosition &gravity() const
        {
            return gravity_;
        }

        GlobalPosition gravity_;
private:
    // the internal method for the initial condition
    // Fully water-saturated saline aquifer: Sw ≈ 1, Sg ≈ 0
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        Scalar pressure_TOP = getParam<Scalar>("Initialization.PressureGWC", 60e5);

        // Compute hydrostatic pressure from top using water density
        FluidState fs;
        fs.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
        fs.setPressure(liquidPhaseIdx, pressure_TOP);
        // Set a pure-water composition to compute liquid density
        fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            if (compIdx != H2OIdx)
                fs.setMoleFraction(liquidPhaseIdx, compIdx, 0.0);

        const Scalar densityW = FluidSystem::density(fs, liquidPhaseIdx);

        // Hydrostatic pressure initialisation (water column)
        priVars[pressureIdx] = pressure_TOP - densityW * gravity()[dimWorld - 1]
                               * (Y_max - globalPos[dimWorld - 1]);

        // TwoPNC: switchIdx = Sn (gas saturation) when both phases present
        // Start fully water-saturated with a numerically safe epsilon gas
        priVars[switchIdx] = 1e-6;

        // Set phase presence state: both phases present (required by primary variable switch)
        priVars.setState(Indices::bothPhases);
        Scalar xgH2O = (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,H2OIdx)* priVars[pressureIdx]) / priVars[pressureIdx];
        priVars[contiCH4EqIdx]  =  1*(1-xgH2O) * priVars[pressureIdx] / (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,CH4Idx) * priVars[pressureIdx]);
        // priVars[contiSO4EqIdx]  =  1.26e-4*(1-xgH2O);
        // Remaining primary variables (dissolved component mole fractions)
        // initialised to zero by the priVars(0.0) constructor - pure brine
        return priVars;
    }
    void initializeOutput()
    {
        Dumux::MetaData::Collector collector;
        Dumux::MetaData::writeJsonFile(collector, name());
    }
        bool enablePressureCutoff_ = true;
        Scalar injPressureMultiplier_ = 1.7;
        Scalar minProdPressure_ = 1e5;
        Scalar pressureInitRef_ = 60e5;
        TimeLoopPtr timeLoop_;
        static constexpr Scalar eps_ = 1e-6;
        bool useNitscheTypeBc_;
        Dimension_vector CELLS_VEC;
        double CELLS_VECX = 0, CELLS_VECY = 0;
        Dimension_vector Upper_Right, Lower_Left;
        Dimension_vector Position0, Position1;
        double Delta_z = 0.0;
        double Delta_x = 0.0;
        double Delta_y = 0.0;
        double X_max = 0.0;
        double X_min = 0;
        double Y_max = 0.0;
        double Y_min = 0.0;
        double Z_max = 0.0;
        double Z_min = 0.0;
        bool enableGravity = false;
        std::string name_;
        // std::string CushionGasType_;
    };

} // end namespace Dumux

#endif
