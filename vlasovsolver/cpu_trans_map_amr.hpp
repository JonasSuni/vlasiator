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
#ifndef CPU_TRANS_MAP_AMR_H
#define CPU_TRANS_MAP_AMR_H

#include <vector>

#include "vec.h"
#include "../common.h"
#include "../spatial_cell.hpp"

struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector< uint > lengthOfPencils; // Lengths of pencils
   std::vector< CellID > ids; // List of cells
   std::vector< Realv > x,y; // x,y - position
   std::vector< bool > periodic;
   std::vector< std::vector<uint> > path; // Path taken through refinement levels

   setOfPencils() {
      
      N = 0;
      sumOfLengths = 0;
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn, bool periodicIn, std::vector<uint> pathIn) {

      N++;
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
      x.push_back(xIn);
      y.push_back(yIn);
      periodic.push_back(periodicIn);
      path.push_back(pathIn);
   }

   std::vector<CellID> getIds(const uint pencilId) const {

      std::vector<CellID> idsOut;
      
      if (pencilId > N) {
         return idsOut;
      }

      CellID ibeg = 0;
      for (uint i = 0; i < pencilId; i++) {
         ibeg += lengthOfPencils[i];
      }
      CellID iend = ibeg + lengthOfPencils[pencilId];
    
      for (uint i = ibeg; i < iend; i++) {
         idsOut.push_back(ids[i]);
      }

      return idsOut;
   }

   // Split one pencil into four pencils covering the same space.
   // dx and dy are the dimensions of the original pencil.
   void split(const uint pencilId, const Realv dx, const Realv dy) {

      auto ids = getIds(pencilId);
      auto path1 = path.at(pencilId);
      auto path2 = path.at(pencilId);
      auto path3 = path.at(pencilId);

      path1.push_back(1);
      path2.push_back(2);
      path3.push_back(3);
      
      x[pencilId] -= 0.25 * dx;
      y[pencilId] += 0.25 * dy;
      path.at(pencilId).push_back(0);     
      
      addPencil(ids, x[pencilId] + 0.25 * dx, y[pencilId] + 0.25 * dy, periodic[pencilId], path1);
      addPencil(ids, x[pencilId] - 0.25 * dx, y[pencilId] - 0.25 * dy, periodic[pencilId], path2);
      addPencil(ids, x[pencilId] + 0.25 * dx, y[pencilId] - 0.25 * dy, periodic[pencilId], path3);
      
   }

};

void compute_spatial_source_cells_for_pencil(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                             setOfPencils pencils,
                                             const uint iPencil,
                                             const uint dimension,
                                             SpatialCell **sourceCells);


void compute_spatial_target_cells_for_pencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                              setOfPencils& pencils,
                                              const uint dimension,
                                              SpatialCell **targetCells);

CellID selectNeighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid,
                      CellID id, int dimension, uint path);

void propagatePencil(Vec* dz, Vec* values, const uint dimension, const uint blockGID,
                     const Realv dt,
                     const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> &vmesh,
                     const uint lengthOfPencil, bool debugflag);

void copy_trans_block_data_amr(
    SpatialCell** source_neighbors,
    const vmesh::GlobalID blockGID,
    int lengthOfPencil,
    Vec* values,
    const unsigned char* const cellid_transpose,
    const uint popID);


setOfPencils buildPencilsWithNeighbors( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid, 
					setOfPencils &pencils, CellID startingId,
					std::vector<CellID> ids, uint dimension, 
					std::vector<uint> path);

void get_seed_ids(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID> &localPropagatedCells,
                  const uint dimension,
                  std::vector<CellID> &seedIds);


bool trans_map_1d_amr(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& localPropagatedCells,
                  const std::vector<CellID>& remoteTargetCells,
                  const uint dimension,
                  const Realv dt,
                  const uint popID);

void update_remote_mapping_contribution(dccrg::Dccrg<spatial_cell::SpatialCell,
                                        dccrg::Cartesian_Geometry>& mpiGrid,
                                        const uint dimension,
                                        int direction,
                                        const uint popID);


#endif
