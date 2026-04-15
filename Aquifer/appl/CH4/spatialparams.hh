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
 * \brief Definition of the spatial parameters for the biogeochemical problem
 *        which uses the isothermal two-phase n-component fully implicit model.
 */



#include <cmath>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCBioGeoTest
 * \brief Definition of the spatial parameters for the biogeochemical problem
 *        which uses the isothermal two-phase n-component fully implicit model.
 */
template<class GridGeometry, class Scalar>
class mixingSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                         mixingSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, mixingSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using Dimension_vector = Dune::FieldVector<Scalar, dimWorld>;
public:
    //! Export the type used for the permeability
    using PermeabilityType = DimWorldMatrix;

    mixingSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry),
    PcKrSwCurve_("SpatialParams.Material")
    {
        referencePorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity", 0.2); //porosity without reactive solid volume fractions
        referencePermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability", 100 * 9.86923e-16);
        temperature_           = getParam<Scalar>("Problem.InitialTemperature");
        DispersionMode         = getParam<int>("Problem.EnableDispersion", 0);
        X_max = this->gridGeometry().bBoxMax()[0];
        X_min = this->gridGeometry().bBoxMin()[0];
        Y_max = this->gridGeometry().bBoxMax()[1];
        Y_min = this->gridGeometry().bBoxMin()[1];
        Z_max = this->gridGeometry().bBoxMax()[2];
        Z_min = this->gridGeometry().bBoxMin()[2];

        // Estimate average cell sizes from grid (approximate for graded grids)
        const auto numElements = this->gridGeometry().gridView().size(0);
        Delta_x = std::cbrt((X_max - X_min) * (Y_max - Y_min) * (Z_max - Z_min) / numElements);
        Delta_y = Delta_x;
        Delta_z = Delta_x;

        if (DispersionMode == 1){
            alphaL_ = Delta_x;
            alphaT_ = Delta_y;
        }
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     */
    template<class ElementSolution>
    Scalar temperature(const Element& element,
                       const SubControlVolume& scv,
                       const ElementSolution& elemSol) const
    { return temperature_; }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        PermeabilityType K(0.0);
        K[0][0] = referencePermeability_;
        K[1][1] = referencePermeability_;
        K[2][2] = referencePermeability_/10;
        return K;
    }



    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return referencePorosity_;
    }

    /*!
     *  \brief Defines the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minimalPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        // const auto eIdx = this->gridGeometry().elementMapper().index(element);
        auto globalPos = element.geometry().center();
        return 1;
        // return M_PI * globalPos[0] * 2;
    }
    
    /*!
     *  \brief Defines the volume fraction of the inert component
     *
     *  \param globalPos The global position in the domain
     *  \param compIdx The index of the inert solid component
     */
    /*!
     *  \brief Defines the reference porosity \f$[-]\f$ distribution.
     *
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return (referencePorosity_); }

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
     * \brief Returns (alphaL, alphaT) for Scheidegger's velocity-aligned dispersion.
     *
     * Scheidegger's model builds the full tensor from velocity direction:
     *   D = (alphaL - alphaT) * v*v^T / |v|  +  alphaT * |v| * I
     * This automatically aligns longitudinal dispersion with flow (radial on
     * a cylindrical grid) — no manual coordinate rotation needed.
     */
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    {
        return { alphaL_, alphaT_ };
    }

    /*!
     * \brief Returns a position-dependent dispersion tensor aligned with the
     *        local radial direction for the FullDispersionTensor model.
     *
     * On a cylindrical grid the "longitudinal" direction at each point is
     * along r = (x, y, 0)/|(x,y)|.  The tensor is rotated so that alphaL
     * applies along r and alphaT in the two perpendicular directions (theta, z).
     */
    DimWorldMatrix dispersionTensor(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    {
        DimWorldMatrix D(0.0);
        // Isotropic base: alphaT in every direction
        for (int i = 0; i < dimWorld; ++i)
            D[i][i] = alphaT_;

        // Add (alphaL - alphaT) in the radial direction: D += (aL-aT) * e_r * e_r^T
        const Scalar rr = std::sqrt(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1]);
        if (rr > 1e-20)
        {
            const Scalar ex = globalPos[0] / rr;
            const Scalar ey = globalPos[1] / rr;
            const Scalar diff = alphaL_ - alphaT_;
            D[0][0] += diff * ex * ex;
            D[0][1] += diff * ex * ey;
            D[1][0] += diff * ey * ex;
            D[1][1] += diff * ey * ey;
            // z-component of radial unit vector is 0, so no z contributions
        }
        return D;
    }

private:
    double Delta_z = 0.0;
    double Delta_x = 0.0;
    double Delta_y = 0.0;
    double X_max = 0.0;
    double X_min = 0;
    double Y_max = 0.0;
    double Y_min = 0.0;
    double Z_max = 0.0;
    double Z_min = 0.0;
    Scalar Permeability_;
    Scalar Porosity_;
    Scalar PermeabilityZ_;
    Scalar referencePorosity_;
    PermeabilityKozenyCarman<PermeabilityType> permLaw_;
    Scalar referencePermeability_;
    PcKrSwCurve PcKrSwCurve_;
    // PorosityPrecipitation<Scalar, /*numFluidComponents*/5, /*activeComponents*/5> poroLaw_;
    Scalar temperature_;
    Scalar initialPyrite_,initialBiofilmMG_,initialBiofilmSRB_;
    Scalar alphaL_ = 0.0;
    Scalar alphaT_ = 0.0;
    int DispersionMode;
};

} // end namespace Dumux
