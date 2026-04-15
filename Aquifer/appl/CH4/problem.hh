// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPNCBioGeoTest
 * \brief Problem for modeling biogeochemical reactions during Underground Hydrogen Storage
 */

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <numeric>
#include <dumux/common/metadata.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/evalgradients.hh>
#include <dumux/io/json/json.hpp>

namespace Dumux
{

    /*!
     * \ingroup TwoPNCBioGeoTest
     * \brief Problem for modeling biogeochemical reactions during Underground Hydrogen Storage
     *
     * To run the simulation execute the following line in shell:
     * <tt>appl_2pncgeobio_biogeochemical_tpfa</tt>
     */
    template <class TypeTag>
    class mixingProblem : public PorousMediumFlowProblem<TypeTag>
    {
        using ParentType = PorousMediumFlowProblem<TypeTag>;
        using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
        using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using Indices = typename ModelTraits::Indices;
        using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;

        // primary variable indices
        enum
        {
            pressureIdx = Indices::pressureIdx,
            switchIdx = Indices::switchIdx
        };

        // phase presence
        enum
        {
            wPhaseOnly = Indices::firstPhaseOnly,
            nPhaseOnly = Indices::secondPhaseOnly,
            bothPhases = Indices::bothPhases
        };

        // equation indices
        enum
        {
            numComponents = FluidSystem::numComponents,

            contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
            contiN2EqIdx = Indices::conti0EqIdx + FluidSystem::N2Idx,
            contiCO2EqIdx = Indices::conti0EqIdx + FluidSystem::CO2Idx,
            contiH2EqIdx = Indices::conti0EqIdx + FluidSystem::H2Idx,
            contiCH4EqIdx = Indices::conti0EqIdx + FluidSystem::CH4Idx,
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


        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
        using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
        using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
        using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
        using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
        using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
        using Element = typename GridView::template Codim<0>::Entity;
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
        using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
        using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
        using SubControlVolume = typename FVElementGeometry::SubControlVolume;
        using FluidState = GetPropType<TypeTag, Properties::FluidState>;
        using SolidState = GetPropType<TypeTag, Properties::SolidState>;
        using ElementIterator = typename GridView::template Codim<0>::Iterator;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
        using ElementGeometry = typename GridGeometry::LocalView;
        using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

        using GridVariablesPtr = std::shared_ptr<GridVariables>;
        GridVariablesPtr gridVariables_;
        enum
        {
            dimWorld = GridView::dimensionworld
        };
        using VelocityOutputType = Dumux::VelocityOutput<GridVariables>;
        using VelocityOutput = VelocityOutputType;
        using Dimension_vector = Dune::FieldVector<Scalar, dimWorld>;
        using Extrusion = Extrusion_t<GridGeometry>;
        //! Property that defines whether mole or mass fractions are used
        static constexpr bool useMoles = ModelTraits::useMoles();
        using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, 1>;
    public:
        mixingProblem(std::shared_ptr<const GridGeometry> gridGeometry)
            : ParentType(gridGeometry), 
            gridVariables_(nullptr)
        {
            if (useMoles)
                std::cout << "problem uses mole-fractions" << std::endl;
            else
                std::cout << "problem uses mass-fractions" << std::endl;

            if (getParam<bool>("Problem.EnableGravity", true))
                gravity_ = {0, 0, -9.81};
            name_ = getParam<std::string>("Problem.Name");
            enableGravity = getParam<bool>("Problem.EnableGravity", true);
            permeability_.resize(gridGeometry->numDofs());
            initializeOutput();

            // Read grid bounding box from actual grid geometry
            X_min = this->gridGeometry().bBoxMin()[0];
            X_max = this->gridGeometry().bBoxMax()[0];
            Y_min = this->gridGeometry().bBoxMin()[1];
            Y_max = this->gridGeometry().bBoxMax()[1];
            Z_min = this->gridGeometry().bBoxMin()[2];
            Z_max = this->gridGeometry().bBoxMax()[2];

            // For cylindrical grid: outer radius from bounding box diagonal in x-y plane
            outerRadius_ = std::max({std::abs(X_max), std::abs(X_min), std::abs(Y_max), std::abs(Y_min)});

            // Well parameters for 3D central well source/sink
            wellRadius_ = getParam<Scalar>("Well.Radius", 0.5);
            wellPerfTop_ = getParam<Scalar>("Well.PerfTop", Z_max);
            wellPerfBottom_ = getParam<Scalar>("Well.PerfBottom", Z_min);
            wellCenterX_ = getParam<Scalar>("Well.CenterX", 0.0);
            wellCenterY_ = getParam<Scalar>("Well.CenterY", 0.0);

            // Precompute total well support volume (sum of cell volumes within well radius)
            wellSupportVolume_ = 0.0;
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                auto fvGeometry = localView(this->gridGeometry());
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& pos = scv.center();
                    const Scalar r = std::sqrt((pos[0] - wellCenterX_)*(pos[0] - wellCenterX_)
                                             + (pos[1] - wellCenterY_)*(pos[1] - wellCenterY_));
                    if (r <= wellRadius_ && pos[2] >= wellPerfBottom_ && pos[2] <= wellPerfTop_)
                        wellSupportVolume_ += scv.volume();
                }
            }
            std::cout << "3D Central Well: radius=" << wellRadius_
                      << " m, supportVolume=" << wellSupportVolume_
                      << " m^3, perf=[" << wellPerfBottom_ << "," << wellPerfTop_ << "]" << std::endl;
             
        }
        const GridVariables& getGridVariables() const {
                return *gridVariables_;
        }
        // double get_time_check() const {return time_check;}
        // void set_time_check(Scalar newValue) const {time_check = newValue;}
        // void LoadGridVariables(GridVariables &gridVariables) { gridvariables_ = gridVariables; } // Public setter method for time_step
        /*!
         * \brief Returns the problem name
         *
         * This is used as a prefix for files generated by the simulation.
         */
        const std::string &name() const
        {
            return name_;
        }

        template<class ElementSolution>
        Scalar extrusionFactor_sub(const Element &element,
                            const SubControlVolume &scv,
                            const ElementSolution& elemSol) const
        {
            return 1; // No extrusion in 3D
        }
        /*!
         * \brief Evaluates the boundary conditions for a Neumann
         *        boundary segment in dependency on the current solution.
         *
         * \param element The finite element
         * \param fvGeometry The finite volume geometry of the element
         * \param elemVolVars All volume variables for the element
         * \param elemFluxVarsCache Flux variables caches for all faces in stencil
         * \param scvf The sub-control volume facefalse
         *
         * This method is used for cases, when the Neumann condition depends on the
         * solution and requires some quantities that are specific to the fully-implicit method.
         * The \a values store the mass flux of each phase normal to the boundary.
         * Negative values indicate an inflow.
         */
        NumEqVector neumann(const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const ElementVolumeVariables &elemVolVars,
                            const ElementFluxVariablesCache &elemFluxVarsCache,
                            const SubControlVolumeFace &scvf) const
        {
            // In 3D, well injection/production is handled via source() term.
            // Neumann boundaries are no-flow (zero flux).
            NumEqVector values(0.0);
            return values;
        }
        void postTimeStep(const SolutionVector &curSol,
                          const SolutionVector &pastSol,
                          const GridVariables &gridVariables,
                          const Scalar timeStepSize,
                          std::shared_ptr<const GridGeometry> fvGridGeometry)
        {
            // Simplified 3D postprocessing: compute inventories and material balance
            Scalar t, dt, volume(0.0);
            NumEqVector inventory(0.0), inventoryPastSol(0.0), materialBalanceError(0.0), inventory_error(0.0);
            Scalar averageReservoirPressure(0.0);

            t = this->time();
            dt = this->timeStepSize();

            for (const auto &element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
            {
                auto fvGeometry = localView(*fvGridGeometry);
                fvGeometry.bind(element);

                auto elemVolVars = localView(gridVariables.curGridVolVars());
                elemVolVars.bindElement(element, fvGeometry, curSol);
                auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
                prevElemVolVars.bindElement(element, fvGeometry, pastSol);

                const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
                const auto elemSolPast = elementSolution(element, pastSol, this->gridGeometry());

                for (auto &&scv : scvs(fvGeometry))
                {
                    const FluidState &fs = elemVolVars[scv].fluidState();
                    const FluidState &fsOld = prevElemVolVars[scv].fluidState();
                    Scalar vol = scv.volume();
                    VolumeVariables volVars, volVarsPast;
                    volVars.update(elemSol, *this, element, scv);
                    volVarsPast.update(elemSolPast, *this, element, scv);

                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        inventory[compIdx] += volVars.porosity() * vol * (
                            fs.saturation(gasPhaseIdx) * fs.molarDensity(gasPhaseIdx) * fs.moleFraction(gasPhaseIdx, compIdx)
                          + fs.saturation(liquidPhaseIdx) * fs.molarDensity(liquidPhaseIdx) * fs.moleFraction(liquidPhaseIdx, compIdx));
                        inventoryPastSol[compIdx] += volVarsPast.porosity() * vol * (
                            fsOld.saturation(gasPhaseIdx) * fsOld.molarDensity(gasPhaseIdx) * fsOld.moleFraction(gasPhaseIdx, compIdx)
                          + fsOld.saturation(liquidPhaseIdx) * fsOld.molarDensity(liquidPhaseIdx) * fsOld.moleFraction(liquidPhaseIdx, compIdx));
                    }
                    averageReservoirPressure += fs.pressure(gasPhaseIdx) * vol;
                    volume += vol;
                }
            }
            const auto &comm = this->gridGeometry().gridView().comm();

            inventory_error = inventoryPastSol - inventory;
            std::cout << "Inventory change (comp 2): " << inventory_error[2] << std::endl;
            std::cout << "---------------------------------------------info end------------------------------------------------------------" << std::endl;

            Dumux::MetaData::Collector collector;
            if (Dumux::MetaData::jsonFileExists(name()))
                Dumux::MetaData::readJsonFile(collector, name());
            collector["time"].push_back(t);
            collector["averageReservoirPressure"].push_back(comm.sum(averageReservoirPressure) / comm.sum(volume));
            for (auto compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                collector["inventory"][FluidSystem::componentName(compIdx)].push_back(comm.sum(inventory[compIdx]));
            }

            Dumux::MetaData::writeJsonFile(collector, name());
        }

        /*!
         * \brief source term: central well injection/production as volumetric source/sink.
         *
         * Well cells are identified by radial distance from (wellCenterX_, wellCenterY_)
         * within wellRadius_ and vertical perforation interval [wellPerfBottom_, wellPerfTop_].
         * Total well rate (mol/s) is distributed uniformly over wellSupportVolume_,
         * so source returns rate per bulk volume [mol/(m^3 s)].
         *
         * Injection: positive source (CH4-focused composition).
         * Production: negative source proportional to local composition.
         */
        NumEqVector source(const Element &element,
                           const FVElementGeometry &fvGeometry,
                           const ElementVolumeVariables &elemVolVars,
                           const SubControlVolume &scv) const
        {
            NumEqVector values(0.0);
            const auto& pos = scv.center();
            const Scalar r = std::sqrt((pos[0] - wellCenterX_)*(pos[0] - wellCenterX_)
                                      + (pos[1] - wellCenterY_)*(pos[1] - wellCenterY_));

            // Only apply source in well cells
            if (r > wellRadius_ || pos[2] < wellPerfBottom_ || pos[2] > wellPerfTop_ || wellSupportVolume_ < 1e-20)
                return values;

            const auto &volVars = elemVolVars[scv];
            Scalar Time = this->time();

            // Operational schedule (same semantics as original code)
            static const Scalar cylcesDev            = getParam<double>("BoundaryConditions.CyclesDev");
            static const Scalar injectionDurationDev = getParam<double>("BoundaryConditions.InjectionDurationDev")*86400;
            static const Scalar idleDurationDev      = getParam<double>("BoundaryConditions.IdleDurationDev")*86400;
            static const Scalar developmentDuration   = cylcesDev*(injectionDurationDev+idleDurationDev);
            static const Scalar idleDurationOp       = getParam<double>("BoundaryConditions.IdleDurationOp")*86400;
            static const Scalar injectionDurationOp  = getParam<double>("BoundaryConditions.InjectionDurationOp")*86400;
            static const Scalar extractionDurationOp = getParam<double>("BoundaryConditions.ExtractionDurationOp")*86400;
            static const Scalar injConc              = getParam<double>("BoundaryConditions.HydrogenInjectionConcentration");
            static const double totalInjRateOp       = getParam<double>("BoundaryConditions.InjectionRateOp", 0.0);
            static const double totalProdRate        = getParam<double>("BoundaryConditions.ProductionRate", 0.0);

            const int cycleNumberOp    = std::floor((Time - developmentDuration)/(injectionDurationOp+extractionDurationOp+idleDurationOp));
            const Scalar localTimeInCycle = Time - developmentDuration - cycleNumberOp*(injectionDurationOp+extractionDurationOp+idleDurationOp);

            // Development phase injection
            if (Time < developmentDuration)
            {
                static const double totalInjRateDev = getParam<double>("BoundaryConditions.InjectionRateDev");
                const int cycleNumber = std::floor(Time/(injectionDurationDev+idleDurationDev));
                const Scalar localTime = Time - cycleNumber*(injectionDurationDev+idleDurationDev);
                if (localTime <= injectionDurationDev)
                {
                    // Inject: positive source = injection into reservoir
                    NumEqVector injComp(0.0);
                    injComp[H2Idx] = injConc;
                    injComp[CO2Idx] = 1.0 - injConc;
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        values[compIdx] = injComp[compIdx] * std::abs(totalInjRateDev) / wellSupportVolume_;
                }
            }
            // Operational phase
            else
            {
                if (localTimeInCycle < injectionDurationOp)
                {
                    // Injection phase: positive source
                    NumEqVector injComp(0.0);
                    injComp[H2Idx] = injConc;
                    injComp[CO2Idx] = 1.0 - injConc;
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        values[compIdx] = injComp[compIdx] * std::abs(totalInjRateOp) / wellSupportVolume_;
                }
                else if (localTimeInCycle > injectionDurationOp + idleDurationOp)
                {
                    // Production phase: negative source proportional to local composition
                    Scalar mobGas = volVars.mobility(gasPhaseIdx);
                    Scalar mobLiq = volVars.mobility(liquidPhaseIdx);
                    Scalar totalMob = mobGas + mobLiq;
                    if (totalMob > 1e-20)
                    {
                        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        {
                            Scalar xGas = volVars.moleFraction(gasPhaseIdx, compIdx);
                            Scalar xLiq = volVars.moleFraction(liquidPhaseIdx, compIdx);
                            // Negative = removal from reservoir
                            values[compIdx] = -std::abs(totalProdRate) / wellSupportVolume_
                                             * (xLiq*mobLiq/totalMob + xGas*mobGas/totalMob);
                        }
                    }
                }
                // else: idle phase, values remain zero
            }

            return values;
        }
    
        /*!
         * \brief Specifies which kind of boundary condition should be
         *        used for which equation on a given boundary segment
         *
         * \param globalPos The global position
         */
        BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
        {
            BoundaryTypes bcTypes;
            // Cylindrical grid: Dirichlet on outer radial boundary, Neumann elsewhere
            const Scalar r = std::sqrt(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1]);
            if (r > outerRadius_ - eps_)
                bcTypes.setAllDirichlet();
            else
                bcTypes.setAllNeumann();
            return bcTypes;
        }

        /*!
         * \brief Evaluates the boundary conditions for a Dirichlet boundary segment
         *
         * \param globalPos The global position
         */
        PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
        {
            return initial_(globalPos);
        }

        

        // Function to return the permeability for additional vtk output
        const std::vector<Scalar> &getPermeability()
        {
            return permeability_;
        }


        // Function to update the permeability for additional vtk output
        template <class SolutionVector>
        void updateVtkOutput(const SolutionVector &curSol)
        {
            for (const auto &element : elements(this->gridGeometry().gridView()))
            {
                const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
                auto fvGeometry = localView(this->gridGeometry());
                fvGeometry.bindElement(element);
                for (const auto &scv : scvs(fvGeometry))
                {
                    VolumeVariables volVars;
                    volVars.update(elemSol, *this, element, scv);
                    const auto dofIdxGlobal = scv.dofIndex();
                    permeability_[dofIdxGlobal] = volVars.permeability()[0][0] / 9.86923e-16;
                }
            }
        }

        /*!
         * \brief Evaluates the initial values for a control volume.
         *
         * \param globalPos The global position
         */
        PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
        {
            return initial_(globalPos);
        }

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

        /*!
         * \brief Evaluates the initial values for a control volume.
         *
         * The internal method for the initial condition
         *
         * \param globalPos The global position
         */

        PrimaryVariables initial_(const GlobalPosition &globalPos) const
        {
            PrimaryVariables values;

            
            // hydrostatic pressure
            
            Scalar pressureGasWaterContact_ = getParam<Scalar>("Initialization.PressureGWC", 100e5);
            Scalar gasWaterContact_ = getParam<Scalar>("Initialization.DepthGWC", 0);
            FluidState fs;
            fs.setMoleFraction(gasPhaseIdx, CO2Idx, 1.0); fs.setMoleFraction(gasPhaseIdx, CH4Idx, 0.0); fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 1);
            const auto &fluidMatrixInteraction = this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
            Scalar thresholdPressure_ = fluidMatrixInteraction.pc(1 - fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr());
            fs.setPressure(liquidPhaseIdx, pressureGasWaterContact_); fs.setPressure(gasPhaseIdx, pressureGasWaterContact_+thresholdPressure_);
            // fs.density(liquidPhaseIdx); fs.density(gasPhaseIdx);
            fs.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
            // Scalar densityW = fs.density(fs, liquidPhaseIdx);
            
            values.setState(bothPhases);

            Scalar pn = 0.0;
            Scalar pw = 0.0;
            Scalar Sw, pC;
            Scalar check = 1.0;
            Scalar z = globalPos[dimWorld - 1];
            Scalar densityW = 1000;//1020;
            Scalar densityG = 150;

            // while (check > 1e-12){
            //     Scalar densityW = FluidSystem::density(fs,liquidPhaseIdx); // 1020;
            //     Scalar densityG = FluidSystem::density(fs,gasPhaseIdx) ;
            //     pw = (pressureGasWaterContact_ - (densityW* gravity()[dimWorld - 1] * (gasWaterContact_ - (globalPos)[dimWorld - 1])));
            //     pn = (pressureGasWaterContact_ + thresholdPressure_ - (densityG * gravity()[dimWorld - 1] * (gasWaterContact_ - (globalPos)[dimWorld - 1])));
            //     fs.setPressure(liquidPhaseIdx, pw); fs.setPressure(gasPhaseIdx, pn);
            //     check = std::abs(densityG - FluidSystem::density(fs,gasPhaseIdx));
            // }
            
            // values[pressureIdx] = pw;
            // values[switchIdx] = 0.8;
            // if (getParam<bool>("Problem.EnableGravity", false)){
            //     values[pressureIdx] = pressureGasWaterContact_;
            //     values[switchIdx] = 1 - fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr();
            //     return values;
            // }
            // Water Zone
            if (z <= gasWaterContact_)
            {
                pw = (pressureGasWaterContact_ - (densityW * gravity()[dimWorld - 1] * (gasWaterContact_ - (globalPos)[dimWorld - 1])));
                pn = pw + thresholdPressure_;

                fs.setPressure(liquidPhaseIdx, pw);
                fs.setPressure(gasPhaseIdx, pn);

                values[pressureIdx] = pw;
                values[switchIdx] = fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr();
            }
            else
            {
                // transition zone_______
                // the water pressure is calculated by the hydrostativ gradient of water
                // this is only true in the transition zone, otherwise it will be overwritten
                while (check > 1e-6){
                    pw = (pressureGasWaterContact_ - (densityW * gravity()[dimWorld - 1] * (gasWaterContact_ - (globalPos)[dimWorld - 1])));
                    // pw = pressureGasWaterContact_ ;
                    fs.setPressure(liquidPhaseIdx, pw);

                    // the gas pressure is calculated by hydrostatic gradient of gas
                    //  this is true in the transition zone and the gas zone
                    pn = (pressureGasWaterContact_ + thresholdPressure_ - (densityG * gravity()[dimWorld - 1] * (gasWaterContact_ - (globalPos)[dimWorld - 1])));
                    // pn = pressureGasWaterContact_ + thresholdPressure_ ;
                    fs.setPressure(gasPhaseIdx, pn);
                    check = std::abs(densityG - FluidSystem::density(fs,gasPhaseIdx))+ std::abs(densityW - FluidSystem::density(fs,liquidPhaseIdx));
                    densityW = FluidSystem::density(fs,liquidPhaseIdx); // 1020;
                    densityG = FluidSystem::density(fs,gasPhaseIdx) ;
                    Sw = fluidMatrixInteraction.pcSwCurve().sw((pn - pw));
                }

                // setting the fluid state in the gas zone
                // for brooks and corey the if statement would never become true, therefore 1% is added
                // which means that the saturation is ver close to the residual water saturation
                if (Sw <= fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr() + 0.0001)
                {
                    Sw = fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr() + 0.0001;
                    values[Indices::switchIdx] = 1 - Sw;

                    // Calculating the water pressure by gas pressure minus capillary pressure
                    pC = fluidMatrixInteraction.pcSwCurve().pc(Sw);
                    pw = pn - pC;
                    fs.setPressure(liquidPhaseIdx, pw);
                    values[Indices::pressureIdx] = pw;
                }
                else
                {
                    values[pressureIdx] = pw;
                    values[switchIdx] = 1 - Sw;
                }
            }
            if (!enableGravity){
                values[pressureIdx] = pressureGasWaterContact_;
                values[switchIdx] = 1 - fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr();
            }

            Scalar xgH2O = (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,H2OIdx)* pw) / pn;
            // values[contiCO2EqIdx]  =  0.5*(1-xgH2O) * pn / (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,CO2Idx) * pw);
            // values[contiCH4EqIdx]  =  0.5*(1-xgH2O) * pn / (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,CH4Idx) * pw);
            // values[contiH2EqIdx]  =  0.5*(1-xgH2O) * pn / (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,H2Idx) * pw);
            // values[contiN2EqIdx]  =  0.1198*(1-xgH2O) * pn / (FluidSystem::fugacityCoefficient(fs,liquidPhaseIdx,N2Idx) * pw);

            return values;
        }

        void initializeOutput()
        {
            Dumux::MetaData::Collector collector;
            Dumux::MetaData::writeJsonFile(collector, name());
        }

        TimeLoopPtr timeLoop_;
        Scalar temperature_;
        DimWorldMatrix preMobility_;
        Scalar mobilityPreLiq;
        Scalar mobilityPreGas;

        Scalar mobilityPostLiq;
        Scalar mobilityPostGas;

        std::vector<double> timeout;
        std::vector<double> bhpout;
        std::vector<double> averageReservoirPressureout;
        Scalar t_old = 0;
        std::vector<double> checkpoints_;
        std::vector<NumEqVector> convertedGeoout, convertedMGout, convertedSRout;
        std::vector<NumEqVector> inventoryout, materialBalanceErrorout;
        std::vector<std::string> componentsOut;
        Scalar initalPyriteContent;
        bool enableMicrobes, enableGeochemistry, enableSulfateReduction;
        
        double X_max = 0.0;
        double X_min = 0.0;
        double Y_max = 0.0;
        double Y_min = 0.0;
        double Z_max = 0.0;
        double Z_min = 0.0;
        static constexpr Scalar eps_ = 1e-6;
        bool enableGravity = false;

        Scalar outerRadius_ = 0.0;

        // 3D central well parameters
        Scalar wellRadius_ = 0.5;
        Scalar wellPerfTop_ = 0.0;
        Scalar wellPerfBottom_ = 0.0;
        Scalar wellCenterX_ = 0.0;
        Scalar wellCenterY_ = 0.0;
        Scalar wellSupportVolume_ = 0.0;

        // the permeability for output
        std::vector<Scalar> permeability_;

        std::string name_;

    };

} // end namespace Dumux

// #endif