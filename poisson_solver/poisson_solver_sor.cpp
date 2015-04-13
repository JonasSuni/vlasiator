/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver_sor.cpp
 * Author: sandroos
 *
 * Created on January 15, 2015, 12:45 PM
 */

#include <cstdlib>
#include <iostream>
#include <omp.h>

#include "../logger.h"
#include "../grid.h"

#include "poisson_solver_sor.h"

using namespace std;

extern Logger logFile;

namespace poisson {

   static const int RED   = 0;
   static const int BLACK = 1;

   vector<CellCache3D> innerCellPointersRED;
   vector<CellCache3D> bndryCellPointersRED;
   vector<CellCache3D> innerCellPointersBLACK;
   vector<CellCache3D> bndryCellPointersBLACK;

   PoissonSolver* makeSOR() {
      return new PoissonSolverSOR();
   }

   PoissonSolverSOR::PoissonSolverSOR(): PoissonSolver() { }

   PoissonSolverSOR::~PoissonSolverSOR() { }

   bool PoissonSolverSOR::initialize() {
      bool success = true;
      return success;
   }

   bool PoissonSolverSOR::finalize() {
      bool success = true;
      return success;
   }

   void PoissonSolverSOR::evaluate(std::vector<poisson::CellCache3D>& cellPointers,const int& cellColor) {
      
      const Real weight = 1.5;

      #pragma omp for
      for (size_t c=0; c<cellPointers.size(); ++c) {
         Real DX2     = cellPointers[c][0][CellParams::DX]; DX2 *= DX2;
         Real DY2     = cellPointers[c][0][CellParams::DY]; DY2 *= DY2;
         Real DZ2     = cellPointers[c][0][CellParams::DZ]; DZ2 *= DZ2;
         Real phi_111 = cellPointers[c][0][CellParams::PHI];
         Real rho_q   = cellPointers[c][0][CellParams::RHOQ_TOT];
         
         Real phi_011 = cellPointers[c][1][CellParams::PHI];
         Real phi_211 = cellPointers[c][2][CellParams::PHI];
         Real phi_101 = cellPointers[c][3][CellParams::PHI];
         Real phi_121 = cellPointers[c][4][CellParams::PHI];
         Real phi_110 = cellPointers[c][5][CellParams::PHI];
         Real phi_112 = cellPointers[c][6][CellParams::PHI];

         Real factor = 2*(1/DX2 + 1/DY2 + 1/DZ2);
         Real rhs = ((phi_011+phi_211)/DX2 + (phi_101+phi_121)/DY2 + (phi_110+phi_112)/DZ2 + rho_q)/factor;
         Real correction = rhs - phi_111;
         cellPointers[c][0][CellParams::PHI] = phi_111 + weight*correction;
      }      
   }

   void cachePointers(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const std::vector<CellID>& cells,std::vector<poisson::CellCache3D>& redCache,
                      std::vector<poisson::CellCache3D>& blackCache) {
      redCache.clear();
      blackCache.clear();
      
      for (size_t c=0; c<cells.size(); ++c) {
         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cells[c]);

         if ((indices[0] + indices[1]%2 + indices[2]%2) % 2 == RED) {
            CellCache3D cache;

            // Cells on domain boundaries are not iterated
            if (mpiGrid[cells[c]]->sysBoundaryFlag != 1) continue;
            
            // Fetch pointers to this cell's (cell) parameters array, 
            // and pointers to +/- xyz face neighbors' arrays
            cache[0] = mpiGrid[cells[c]]->parameters;
            
            indices[0] -= 1; cache[1] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] += 2; cache[2] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] -= 1;
            
            indices[1] -= 1; cache[3] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] += 2; cache[4] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] -= 1;
            
            indices[2] -= 1; cache[5] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] += 2; cache[6] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] -= 1;
            
            redCache.push_back(cache);
         } else {
            CellCache3D cache;
            
            // Cells on domain boundaries are not iterated
            if (mpiGrid[cells[c]]->sysBoundaryFlag != 1) continue;
            
            // Fetch pointers to this cell's (cell) parameters array,
            // and pointers to +/- xyz face neighbors' arrays
            cache[0] = mpiGrid[cells[c]]->parameters;
            
            indices[0] -= 1; cache[1] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] += 2; cache[2] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] -= 1;
            
            indices[1] -= 1; cache[3] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] += 2; cache[4] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] -= 1;
            
            indices[2] -= 1; cache[5] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] += 2; cache[6] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] -= 1;
            
            blackCache.push_back(cache);
         }
      }
   }

   bool PoissonSolverSOR::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;

      // If mesh partitioning has changed, recalculate pointer caches
      if (Parameters::meshRepartitioned == true) {
         phiprof::start("Pointer Caching");
         cachePointers(mpiGrid,mpiGrid.get_local_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID),bndryCellPointersRED,bndryCellPointersBLACK);
         cachePointers(mpiGrid,mpiGrid.get_local_cells_not_on_process_boundary(POISSON_NEIGHBORHOOD_ID),innerCellPointersRED,innerCellPointersBLACK);
         phiprof::stop("Pointer Caching");
      }

      // Update charge density
      phiprof::start("MPI (RHOQ)");
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
      phiprof::stop("MPI (RHOQ)");

      do {
         int iterations = 0;
         const int N_iterations = 10;

	 #pragma omp parallel
	   {
	      const int tid = omp_get_thread_num();
	      
	      // Iterate the potential N_iterations times and then
	      // check if the error is less than the required value
	      for (int N=0; N<N_iterations; ++N) {
		 // Make a copy of the potential if we are going 
		 // to evaluate the solution error
		 if (N == N_iterations-1) {
		    if (tid == 0) phiprof::start("Copy Old Potential");
		    #pragma omp for
		    for (size_t c=0; c<Poisson::localCellParams.size(); ++c) {
		       Poisson::localCellParams[c][CellParams::PHI_TMP] = Poisson::localCellParams[c][CellParams::PHI];
		    }
		    if (tid == 0) phiprof::stop("Copy Old Potential");
		 }

		 // Solve red cells first, the black cells
		 if (solve(mpiGrid,RED  ) == false) success = false;
		 if (solve(mpiGrid,BLACK) == false) success = false;
	      }
	   }

         // Evaluate the error in potential solution and reiterate if necessary
         iterations += N_iterations;
         const Real relPotentialChange = error(mpiGrid);
         if (relPotentialChange <= Poisson::minRelativePotentialChange) break;
         if (iterations >= Poisson::maxIterations) break;
      } while (true);

      return success;
   }

   bool PoissonSolverSOR::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const int& oddness) {
      // NOTE: This function is entered by all threads in OpenMP run,
      // so everything must be thread-safe!

      bool success = true;
      const int tid = omp_get_thread_num();

      // Compute new potential on process boundary cells
      if (tid == 0) phiprof::start("Evaluate potential");
      if (oddness == RED) evaluate(bndryCellPointersRED,oddness);
      else                evaluate(bndryCellPointersBLACK,oddness);
      if (tid == 0) {
	 phiprof::stop("Evaluate potential");
	 
	 // Exchange new potential values on process boundaries
	 phiprof::start("MPI (start copy)");
	 mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
	 phiprof::stop("MPI (start copy)");
	 
	 phiprof::start("Evaluate potential");
      }

      // Compute new potential on inner cells
      if (oddness == RED) evaluate(innerCellPointersRED,oddness);
      else                evaluate(innerCellPointersBLACK,oddness);
      
      // Wait for MPI transfers to complete
      if (tid == 0) {
	 phiprof::stop("Evaluate potential");
	 phiprof::start("MPI (wait copy)");
	 mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
	 phiprof::stop("MPI (wait copy)");
      }

      return success;
   }
   
} // namespace poisson
