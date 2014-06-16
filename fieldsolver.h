/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#ifndef FIELDSOLVER_H
#define FIELDSOLVER_H

#include "definitions.h"
#include "common.h"
#include "spatial_cell.hpp"
#include "sysboundary/sysboundary.h"
using namespace spatial_cell;

#include <dccrg.hpp>


/*
namespace fieldsolver {
   
   enum RecVars {drhodx,drhody,drhodz,
	dBxdy,dBxdz,dBydx,dBydz,dBzdx,dBzdy,
	dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz
   };
   
} // namespace fieldsolver
*/

void calculateVolumeAveragedFields(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
bool finalizeFieldPropagator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
bool initializeFieldPropagator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                               SysBoundary& sysBoundaries);
bool initializeFieldPropagatorAfterRebalance(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
bool propagateFields(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                     SysBoundary& sysBoundaries,
                     creal& dt);
void calculateEdgeElectricFieldX(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);
void calculateEdgeElectricFieldY(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);
void calculateEdgeElectricFieldZ(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);

CellID getNeighbourID(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k
);
Real limiter(creal& left,creal& cent,creal& rght);
Real divideIfNonZero(creal rhoV, creal rho);
#endif
