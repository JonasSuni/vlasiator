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

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "GEMchallenge.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   GEMchallenge::GEMchallenge(): TriAxisSearch() { }
   GEMchallenge::~GEMchallenge() { }

   bool GEMchallenge::initialize(void) {return Project::initialize();}

   void GEMchallenge::addParameters(){
      typedef Readparameters RP;
      RP::add("GEMchallenge.Scale_size", "GEMchallenge sheet scale size (m)", 150000.0);
      RP::add("GEMchallenge.VX0", "Initial Velocity in x-direction", 1e4);
      RP::add("GEMchallenge.BX0", "Magnetic field at infinity (T)", 1);
      RP::add("GEMchallenge.BY0", "Magnetic field at infinity (T)", 0.0);
      RP::add("GEMchallenge.BZ0", "Magnetic field at infinity (T)", 0.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_GEMchallenge.Temperature", "Temperature (K)", 2.0e6);
         RP::add(pop + "_GEMchallenge.rho", "Number density at infinity (m^-3)", 1.0e7);
      }
   }

   void GEMchallenge::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("GEMchallenge.Scale_size", this->SCA_LAMBDA);
      RP::get("GEMchallenge.VX0", this->VX0);
      RP::get("GEMchallenge.BX0", this->BX0);
      RP::get("GEMchallenge.BY0", this->BY0);
      RP::get("GEMchallenge.BZ0", this->BZ0);


      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         GEMchallengeSpeciesParameters sP;

         RP::get(pop + "_GEMchallenge.Temperature", sP.TEMPERATURE);
         RP::get(pop + "_GEMchallenge.rho", sP.DENSITY);

         speciesParams.push_back(sP);
      }
   }

   Realf GEMchallenge::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const GEMchallengeSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
      
      // creal kz = 8 * M_PI / (Parameters::zmax - Parameters::zmin);

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY / pow(cosh(z / (this->SCA_LAMBDA)), 2.0) + sP.DENSITY * 0.2;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      // creal rhofac = (this->BX0*this->BX0 + this->BY0*this->BY0 + this->BZ0*this->BZ0) / 2.0 / physicalconstants::MU_0  / physicalconstants::K_B / initT;

      // initRho += rhofac / pow(cosh(z / (this->SCA_LAMBDA)), 2.0);

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
               const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

   /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf GEMchallenge::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const GEMchallengeSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY / pow(cosh(z / (this->SCA_LAMBDA)), 2.0) + sP.DENSITY * 0.2;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      // creal rhofac = (this->BX0*this->BX0 + this->BY0*this->BY0 + this->BZ0*this->BZ0) / 2.0 / physicalconstants::MU_0  / physicalconstants::K_B / initT;

      // initRho += rhofac / pow(cosh(z / (this->SCA_LAMBDA)), 2.0);
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }

   void GEMchallenge::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   vector<std::array<Real, 3>> GEMchallenge::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> V0;
      std::array<Real, 3> v = {{0.0, 0.0, 0.0 }};
      V0.push_back(v);
      return V0;
   }

   void GEMchallenge::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);

      Real Bx_island, By_island, Bz_island;
      creal Lx = Parameters::xmax - Parameters::xmin;
      creal Ly = Parameters::ymax - Parameters::ymin;
      creal Lz = Parameters::zmax - Parameters::zmin;

      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

                  Bx_island = -M_PI * this->BX0 * 0.1 * Lz * cos(2.0 * M_PI * (xyz[0] + 0.5 * perBGrid.DX) / Lx) * sin(M_PI * (xyz[2] + 0.5 * perBGrid.DZ) / Lz);
                  Bz_island = 2.0 * M_PI * Lx * this->BX0 * 0.1 * sin(2.0 * M_PI * (xyz[0] + 0.5 * perBGrid.DX) / Lx) * cos(M_PI * (xyz[2] + 0.5 * perBGrid.DZ) / Lz);

                  cell->at(fsgrids::bfield::PERBX) = this->BX0 * tanh((xyz[1] + 0.5 * perBGrid.DY) / this->SCA_LAMBDA) + Bx_island;
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = 0.0 + Bz_island;
               }
            }
         }
      }
   }

} // namespace projects
