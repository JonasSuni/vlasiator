/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "fs_common.h"
#include "ldz_hall.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

// Template function for Hall term calculations
template<typename REAL>
inline REAL calculateHallTermComponent(
    const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC,
    creal bgField1, creal bgField2,
    creal dx, creal dy, creal dz,
    HallComponent component,
    bool isPositive  // Determines sign of coefficients
) {
    // Common calculations for all components
    REAL result = 0.0;
    
    switch(component) {
        case HallComponent::X:
            // Consolidated JXBX calculations
            result = isPositive ? 
                JXBX_000_100(pC, bgField1, bgField2, dx, dy, dz) :
                JXBX_010_110(pC, bgField1, bgField2, dx, dy, dz);
            break;
        case HallComponent::Y:
            // Consolidated JXBY calculations
            result = isPositive ?
                JXBY_000_010(pC, bgField1, bgField2, dx, dy, dz) :
                JXBY_100_110(pC, bgField1, bgField2, dx, dy, dz);
            break;
        case HallComponent::Z:
            // Consolidated JXBZ calculations
            result = isPositive ?
                JXBZ_000_001(pC, bgField1, bgField2, dx, dy, dz) :
                JXBZ_100_101(pC, bgField1, bgField2, dx, dy, dz);
            break;
    }
    return result;
}

// Optimized edge Hall term calculation
template<HallComponent C>
void calculateEdgeHallTermComponents(
    // ...existing parameters...
) {
    if (Parameters::ohmHallTerm == 0) {
        cerr << __FILE__ << ":" << __LINE__ << " Invalid Hall term configuration" << endl;
        return;
    }

    const Real dx = technicalGrid.DX;
    const Real dy = technicalGrid.DY;
    const Real dz = technicalGrid.DZ;

    // Get appropriate background fields based on component
    Real bgField1, bgField2;
    if constexpr (C == HallComponent::X) {
        bgField1 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY);
        bgField2 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ);
    } else if constexpr (C == HallComponent::Y) {
        bgField1 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX);
        bgField2 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ);
    } else {
        bgField1 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX);
        bgField2 = BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY);
    }

    if (Parameters::ohmHallTerm == 1) {
        // Simple Hall term calculation
        // ...existing simple calculation code...
    } else {
        // Advanced Hall term calculation using template
        std::array<Real,4> hallRhoq;
        
        // Calculate charge densities for all corners
        for(int corner = 0; corner < 4; corner++) {
            const int di = (corner & 1) ? 1 : -1;
            const int dj = (corner & 2) ? 1 : -1;
            
            hallRhoq[corner] = FOURTH * (
                momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) +
                momentsGrid.get(i+di,j,k)->at(fsgrids::moments::RHOQ) +
                momentsGrid.get(i,j+dj,k)->at(fsgrids::moments::RHOQ) +
                momentsGrid.get(i+di,j+dj,k)->at(fsgrids::moments::RHOQ)
            );
            
            hallRhoq[corner] = std::max(hallRhoq[corner], Parameters::hallMinimumRhoq);
        }

        // Calculate Hall terms for each corner
        const auto hallTerms = calculateHallTermForAllCorners<C>(
            perturbedCoefficients, bgField1, bgField2, 
            dx, dy, dz, hallRhoq
        );

        // Store results in appropriate components
        storeHallTermResults<C>(EHallGrid, i, j, k, hallTerms);
    }
}

// Main Hall term calculation
void calculateHallTerm(
    // ...existing parameters...
) {
    // ...existing validation code...

    std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients;
    reconstructionCoefficients(perBGrid, dPerBGrid, perturbedCoefficients, 
                             i, j, k, 3);

    if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && 
        (cellSysBoundaryLayer != 1)) {
        // Handle boundary conditions
        for(int component = 0; component < 3; component++) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
                ->fieldSolverBoundaryCondHallElectricField(
                    EHallGrid, i, j, k, component);
        }
    } else {
        // Calculate components
        calculateEdgeHallTermComponents<HallComponent::X>(
            perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid,
            BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
        calculateEdgeHallTermComponents<HallComponent::Y>(
            perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid,
            BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
        calculateEdgeHallTermComponents<HallComponent::Z>(
            perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid,
            BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
    }
}

// ...existing calculateHallTermSimple implementation...
