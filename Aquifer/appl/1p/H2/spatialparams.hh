// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for the 2pnc problems.
 */

#ifndef DUMUX_2PNC_TEST_SPATIAL_PARAMS_HH
#define DUMUX_2PNC_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for the 2pnc test problems.
 */
template<class GridGeometry, class Scalar>
class TwoPNCTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                             TwoPNCTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           TwoPNCTestSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using Dimension_vector = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    TwoPNCTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry),
      PcKrSwCurve_("SpatialParams.Material")
    {
        referencePorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity", 0.2);
        referencePermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability", 100 * 9.86923e-16);
        temperature_           = getParam<Scalar>("Problem.InitialTemperature");
        DispersionMode         = getParam<int>("Problem.DispersionMode", 0);
        CELLS_VEC = getParam<Dimension_vector>("Grid.Cells");
        X_max = this->gridGeometry().bBoxMax()[0];
        X_min = this->gridGeometry().bBoxMin()[0];
        Y_max = this->gridGeometry().bBoxMax()[1];
        Y_min = this->gridGeometry().bBoxMin()[1];

        Delta_x = (X_max - X_min) / CELLS_VEC[0];
        Delta_y = (Y_max - Y_min) / CELLS_VEC[1];
        if (DispersionMode == 1){
            alphaL_ = Delta_x;
            alphaT_ = Delta_y;
        }

        for (int i = 0; i < dimWorld; i++){
            if(i == 0)
                dispersionTensor_[i][i] = alphaL_;
            else
                dispersionTensor_[i][i] = alphaT_;
        }
    }

    template<class ElementSolution>
    Scalar temperature(const Element& element,
                       const SubControlVolume& scv,
                       const ElementSolution& elemSol) const
    { return temperature_; }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        PermeabilityType K(0.0);
        K[0][0] = referencePermeability_;
        K[1][1] = referencePermeability_/10;
        return K;
    }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return referencePorosity_;
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return referencePorosity_; }

    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }

    Scalar minimalPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(PcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    /*!
     * \brief Defines the dispersion tensor \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    const DimWorldMatrix &dispersionTensor(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    {
        return dispersionTensor_;
    }

private:
    Dimension_vector CELLS_VEC;
    double Delta_x = 0.0;
    double Delta_y = 0.0;
    double X_max = 0.0;
    double X_min = 0.0;
    double Y_max = 0.0;
    double Y_min = 0.0;
    Scalar referencePorosity_;
    Scalar referencePermeability_;
    PcKrSwCurve PcKrSwCurve_;
    Scalar temperature_;
    Scalar alphaL_ = 0.0;
    Scalar alphaT_ = 0.0;
    DimWorldMatrix dispersionTensor_;
    int DispersionMode;
};

} // end namespace Dumux

#endif
