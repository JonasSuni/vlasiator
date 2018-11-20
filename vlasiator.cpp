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
#include <vector>
#include <sstream>
#include <ctime>
#include <omp.h>

#ifdef _OPENMP
   #include <omp.h>
#endif

#include <fsgrid.hpp>

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"
#include "readparameters.h"
#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"

#include "fieldsolver/fs_common.h"
#include "poisson_solver/poisson_solver.h"
#include "projects/project.h"
#include "grid.h"
#include "iowrite.h"
#include "ioread.h"

#include "object_wrapper.h"
#include "fieldsolver/gridGlue.hpp"

#ifdef CATCH_FPE
#include <fenv.h>
#include <signal.h>
/*! Function used to abort the program upon detecting a floating point exception. Which exceptions are caught is defined using the function feenableexcept.
 */
void fpehandler(int sig_num)
{
   signal(SIGFPE, fpehandler);
   printf("SIGFPE: floating point exception occured, exiting.\n");
   abort();
}
#endif

#include "phiprof.hpp"

Logger logFile, diagnostic;
static dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> mpiGrid;

using namespace std;
using namespace phiprof;

int globalflags::bailingOut = 0;
bool globalflags::writeRestart = 0;
bool globalflags::balanceLoad = 0;

ObjectWrapper objectWrapper;

void addTimedBarrier(string name){
#ifdef NDEBUG
//let's not do  a barrier
   return; 
#endif
   int bt=phiprof::initializeTimer(name,"Barriers","MPI");
   phiprof::start(bt);
   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::stop(bt);
}

bool computeNewTimeStep(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,Real &newDt, bool &isChanged) {
   
   phiprof::start("compute-timestep");
   //compute maximum time-step, this cannot be done at the first
   //step as the solvers compute the limits for each cell

   isChanged=false;

   const vector<CellID>& cells = getLocalCells();
   /* Arrays for storing local (per process) and global max dt
      0th position stores ordinary space propagation dt
      1st position stores velocity space propagation dt
      2nd position stores field propagation dt
   */
   Real dtMaxLocal[3];
   Real dtMaxGlobal[3];
   
   dtMaxLocal[0]=numeric_limits<Real>::max();
   dtMaxLocal[1]=numeric_limits<Real>::max();
   dtMaxLocal[2]=numeric_limits<Real>::max();

   for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
      SpatialCell* cell = mpiGrid[*cell_id];
      const Real dx = cell->parameters[CellParams::DX];
      const Real dy = cell->parameters[CellParams::DY];
      const Real dz = cell->parameters[CellParams::DZ];
      
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         const Real* blockParams = blockContainer.getParameters();
         const Real EPS = numeric_limits<Real>::min()*1000;
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            for (unsigned int i=0; i<WID;i+=WID-1) {
                const Real Vx 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VXCRD] 
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVX]
                  + EPS;
                const Real Vy 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VYCRD] 
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVY]
                  + EPS;
                const Real Vz 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VZCRD]
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVZ]
                  + EPS;

                const Real dt_max_cell = min(dx/fabs(Vx),min(dy/fabs(Vy),dz/fabs(Vz)));
                cell->parameters[CellParams::MAXRDT] = min(dt_max_cell,cell->parameters[CellParams::MAXRDT]);
                cell->set_max_r_dt(popID,min(dt_max_cell,cell->get_max_r_dt(popID)));
             }
         }
      }
      
      
      if ( cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
           (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )) {
         //spatial fluxes computed also for boundary cells
         dtMaxLocal[0]=min(dtMaxLocal[0], cell->parameters[CellParams::MAXRDT]);
         dtMaxLocal[2]=min(dtMaxLocal[2], cell->parameters[CellParams::MAXFDT]);
      }

      if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY && cell->parameters[CellParams::MAXVDT] != 0) {
         //Acceleration only done on non sysboundary cells
         dtMaxLocal[1]=min(dtMaxLocal[1], cell->parameters[CellParams::MAXVDT]);
      }
   }
   MPI_Allreduce(&(dtMaxLocal[0]), &(dtMaxGlobal[0]), 3, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
   
   //If any of the solvers are disabled there should be no limits in timespace from it
   if (P::propagateVlasovTranslation == false)
      dtMaxGlobal[0]=numeric_limits<Real>::max();
   if (P::propagateVlasovAcceleration == false)
      dtMaxGlobal[1]=numeric_limits<Real>::max();
   if (P::propagateField == false)
      dtMaxGlobal[2]=numeric_limits<Real>::max();
   
   creal meanVlasovCFL = 0.5*(P::vlasovSolverMaxCFL+ P::vlasovSolverMinCFL);
   creal meanFieldsCFL = 0.5*(P::fieldSolverMaxCFL+ P::fieldSolverMinCFL);
   Real subcycleDt;

   // cout << dtMaxGlobal[0];
   // cout << " ";
   // cout << dtMaxGlobal[1];
   // cout << " ";
   // cout << dtMaxGlobal[2];
   // cout << " ";
   // cout << P::vlasovSolverMaxCFL;
   // cout << " ";
   // cout << P::vlasovSolverMinCFL;
   // cout << " ";
   // cout << P::fieldSolverMaxCFL;
   // cout << " ";
   // cout << P::fieldSolverMinCFL;
   // cout << endl;
   
   
   //reduce dt if it is too high for any of the three propagators, or too low for all propagators
   if(( P::dt > dtMaxGlobal[0] * P::vlasovSolverMaxCFL ||
        P::dt > dtMaxGlobal[1] * P::vlasovSolverMaxCFL * P::maxSlAccelerationSubcycles ||
        P::dt > dtMaxGlobal[2] * P::fieldSolverMaxCFL * P::maxFieldSolverSubcycles ) ||
      ( P::dt < dtMaxGlobal[0] * P::vlasovSolverMinCFL && 
        P::dt < dtMaxGlobal[1] * P::vlasovSolverMinCFL * P::maxSlAccelerationSubcycles &&
        P::dt < dtMaxGlobal[2] * P::fieldSolverMinCFL * P::maxFieldSolverSubcycles )
      ) {

      //new dt computed
      isChanged=true;

      //set new timestep to the lowest one of all interval-midpoints
      const Real half = 0.5;
      newDt = meanVlasovCFL * dtMaxGlobal[0];
      newDt = min(newDt,meanVlasovCFL * dtMaxGlobal[1] * P::maxSlAccelerationSubcycles);
      newDt = min(newDt,meanFieldsCFL * dtMaxGlobal[2] * P::maxFieldSolverSubcycles);
   
      logFile <<"(TIMESTEP) New dt = " << newDt << " computed on step "<<  P::tstep <<" at " <<P::t << 
         "s   Maximum possible dt (not including  vlasovsolver CFL "<< 
         P::vlasovSolverMinCFL <<"-"<<P::vlasovSolverMaxCFL<<
         " or fieldsolver CFL "<< 
         P::fieldSolverMinCFL <<"-"<<P::fieldSolverMaxCFL<<
         ") in {r, v, BE} was " <<
         dtMaxGlobal[0] << " " <<
         dtMaxGlobal[1] << " " <<
         dtMaxGlobal[2] << " " <<
         " Including subcycling { v, BE}  was " <<
         dtMaxGlobal[1] * P::maxSlAccelerationSubcycles << " " <<
         dtMaxGlobal[2] * P::maxFieldSolverSubcycles<< " " <<
         endl << writeVerbose;

      if(P::dynamicTimestep == true) {
         subcycleDt = newDt;
      } else {
         logFile <<"(TIMESTEP) However, fixed timestep in config overrides dt = " << P::dt << endl << writeVerbose;
         subcycleDt = P::dt;
      }
   } else {
      subcycleDt = P::dt;
   }
   
   // Subcycle if field solver dt < global dt (including CFL) (new or old dt hence the hassle with subcycleDt
   if (meanFieldsCFL*dtMaxGlobal[2] < subcycleDt && P::propagateField) {
      P::fieldSolverSubcycles = min(convert<uint>(ceil(subcycleDt / (meanFieldsCFL*dtMaxGlobal[2]))), P::maxFieldSolverSubcycles);
   } else {
      P::fieldSolverSubcycles = 1;
   }
   
   phiprof::stop("compute-timestep");
   return true;
}

ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}

/** Get local cell IDs. This function creates a cached copy of the 
 * cell ID lists to significantly improve performance. The cell ID 
 * cache is recalculated every time the mesh partitioning changes.
 * @return Local cell IDs.*/
const std::vector<CellID>& getLocalCells() {
   return Parameters::localCells;
}

void recalculateLocalCellsCache() {
     {
        vector<CellID> dummy;
        dummy.swap(Parameters::localCells);
     }
   Parameters::localCells = mpiGrid.get_cells();
}

int main(int argn,char* args[]) {
   bool success = true;
   int myRank, doBailout;
   const creal DT_EPSILON=1e-12;
   typedef Parameters P;
   Real newDt;
   bool dtIsChanged;
   
   const bool printLines = false;
   const bool printCells = false;
   const bool printSums  = true;
   
// Init MPI:
   int required=MPI_THREAD_FUNNELED;
   int provided;
   MPI_Init_thread(&argn,&args,required,&provided);
   if (required > provided){
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(myRank==MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed! Got " << provided << ", need "<<required <<endl;
      exit(1);
   }
   
   phiprof::initialize();
   
   double initialWtime =  MPI_Wtime();
   
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm,&myRank);
   SysBoundary sysBoundaries;
   bool isSysBoundaryCondDynamic;
   
   #ifdef CATCH_FPE
   // WARNING FE_INEXACT is too sensitive to be used. See man fenv.
   //feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
   feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
   //feenableexcept(FE_DIVBYZERO|FE_INVALID);
   signal(SIGFPE, fpehandler);
   #endif

   phiprof::start("main");
   phiprof::start("Initialization");
   phiprof::start("Read parameters");
   //init parameter file reader
   Readparameters readparameters(argn,args,MPI_COMM_WORLD);
   P::addParameters();
   getObjectWrapper().addParameters();
   readparameters.parse(); // First pass parsing
   if (P::getParameters() == false) {
      if (myRank == MASTER_RANK) {
         cerr << "(MAIN) ERROR: getParameters failed!" << endl;
      }
      exit(1);
   }

   getObjectWrapper().addPopulationParameters();
   sysBoundaries.addParameters();
   projects::Project::addParameters();
   Project* project = projects::createProject();
   getObjectWrapper().project = project;
   
   readparameters.parse(); // Second pass parsing: specific population parameters
   readparameters.helpMessage(); // Call after last parse, exits after printing help if help requested
   getObjectWrapper().getParameters();
   project->getParameters();
   sysBoundaries.getParameters();
   phiprof::stop("Read parameters");

   // Init parallel logger:
   phiprof::start("open logFile & diagnostic");
   //if restarting we will append to logfiles
   if (logFile.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt",P::isRestart) == false) {
      if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open logfile!" << endl;
      exit(1);
   }
   if (P::diagnosticInterval != 0) {
      if (diagnostic.open(MPI_COMM_WORLD,MASTER_RANK,"diagnostic.txt",P::isRestart) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open diagnostic file!" << endl;
         exit(1);
      }
   }
   {
      int mpiProcs;
      MPI_Comm_size(MPI_COMM_WORLD,&mpiProcs);
      logFile << "(MAIN) Starting simulation with " << mpiProcs << " MPI processes ";
      #ifdef _OPENMP
         logFile << "and " << omp_get_max_threads();
      #else
         logFile << "and 0";
      #endif
      logFile << " OpenMP threads per process" << endl << writeVerbose;      
   }
   phiprof::stop("open logFile & diagnostic");
   
   // Init project
   phiprof::start("Init project");
   if (project->initialize() == false) {
      if(myRank == MASTER_RANK) cerr << "(MAIN): Project did not initialize correctly!" << endl;
      exit(1);
   }
   if (project->initialized() == false) {
      if (myRank == MASTER_RANK) {
         cerr << "(MAIN): Project base class was not initialized!" << endl;
         cerr << "\t Call Project::initialize() in your project's initialize()-function." << endl;
         exit(1);
      }
   }
   phiprof::stop("Init project");
   
   // Add AMR refinement criterias:
   amr_ref_criteria::addRefinementCriteria();

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   // Initialize grid.  After initializeGrid local cells have dist
   // functions, and B fields set. Cells have also been classified for
   // the various sys boundary conditions.  All remote cells have been
   // created. All spatial date computed this far is up to date for
   // FULL_NEIGHBORHOOD. Block lists up to date for
   // VLASOV_SOLVER_NEIGHBORHOOD (but dist function has not been communicated)
   phiprof::start("Init grid");
   //dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> mpiGrid;
   initializeGrid(argn,args,mpiGrid,sysBoundaries,*project);
   isSysBoundaryCondDynamic = sysBoundaries.isDynamic();

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   phiprof::stop("Init grid");
   
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   phiprof::start("Init DROs");
   DataReducer outputReducer, diagnosticReducer;
   initializeDataReducers(&outputReducer, &diagnosticReducer);
   phiprof::stop("Init DROs");  
   
   // Initialize simplified Fieldsolver grids.
   phiprof::start("Init fieldsolver grids");
   const std::array<int,3> fsGridDimensions = {convert<int>(P::xcells_ini) * pow(2,P::amrMaxSpatialRefLevel),
                                               convert<int>(P::ycells_ini) * pow(2,P::amrMaxSpatialRefLevel),
                                               convert<int>(P::zcells_ini) * pow(2,P::amrMaxSpatialRefLevel)};

   std::array<bool,3> periodicity{mpiGrid.topology.is_periodic(0),
                                  mpiGrid.topology.is_periodic(1),
                                  mpiGrid.topology.is_periodic(2)};

   const int fsGridSize = (fsGridDimensions[0] + 4) * (fsGridDimensions[1] + 4) * (fsGridDimensions[2] + 4);
   // adding 1 just to be safe
   const int tagOffset = fsGridSize + 1;
   int tagId = 0;
   
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> perBGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> perBDt2Grid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> EGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> EDt2Grid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> EHallGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> EGradPeGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> momentsGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> momentsDt2Grid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> dPerBGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> dMomentsGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> BgBGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> volGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId++);
   FsGrid< fsgrids::technical, 2> technicalGrid(fsGridDimensions, comm, periodicity, tagOffset * tagId);
   // Set DX,DY and DZ
   // TODO: This is currently just taking the values from cell 1, and assuming them to be
   // constant throughout the simulation.
   perBGrid.DX = perBDt2Grid.DX = EGrid.DX = EDt2Grid.DX = EHallGrid.DX = EGradPeGrid.DX = momentsGrid.DX
      = momentsDt2Grid.DX = dPerBGrid.DX = dMomentsGrid.DX = BgBGrid.DX = volGrid.DX = technicalGrid.DX
      = P::dx_ini * pow(2,-P::amrMaxSpatialRefLevel);
   perBGrid.DY = perBDt2Grid.DY = EGrid.DY = EDt2Grid.DY = EHallGrid.DY = EGradPeGrid.DY = momentsGrid.DY
      = momentsDt2Grid.DY = dPerBGrid.DY = dMomentsGrid.DY = BgBGrid.DY = volGrid.DY = technicalGrid.DY
      = P::dy_ini * pow(2,-P::amrMaxSpatialRefLevel);
   perBGrid.DZ = perBDt2Grid.DZ = EGrid.DZ = EDt2Grid.DZ = EHallGrid.DZ = EGradPeGrid.DZ = momentsGrid.DZ
      = momentsDt2Grid.DZ = dPerBGrid.DZ = dMomentsGrid.DZ = BgBGrid.DZ = volGrid.DZ = technicalGrid.DZ
      = P::dz_ini * pow(2,-P::amrMaxSpatialRefLevel);
   phiprof::stop("Init fieldsolver grids");
   phiprof::start("Initial fsgrid coupling");
   const std::vector<CellID>& cells = getLocalCells();

   // for (auto cell: cells) {
   //    creal dx = mpiGrid[cell]->parameters[CellParams::DX];
   //    creal dy = mpiGrid[cell]->parameters[CellParams::DY];
   //    creal dz = mpiGrid[cell]->parameters[CellParams::DZ];

   //    std::cout << "cell " << cell << " dx,dy,dz = " << dx << ", " << dy << ", " << dz << std::endl;
   // }

   // cout << "Local cells are: ";
   // for(auto id : cells) cout << id << " ";
   // cout << endl;
   
   // Couple FSGrids to mpiGrid
   // TODO: Do we really need to couple *all* of these fields?

   bool debugFsgrid = false;
   
   perBGrid.      setupForGridCoupling(debugFsgrid);   
   perBDt2Grid.   setupForGridCoupling(debugFsgrid);      
   EGrid.         setupForGridCoupling(debugFsgrid);
   EDt2Grid.      setupForGridCoupling(debugFsgrid);
   EHallGrid.     setupForGridCoupling(debugFsgrid);
   EGradPeGrid.   setupForGridCoupling(debugFsgrid);
   momentsGrid.   setupForGridCoupling(debugFsgrid);
   momentsDt2Grid.setupForGridCoupling(debugFsgrid);
   dPerBGrid.     setupForGridCoupling(debugFsgrid);
   dMomentsGrid.  setupForGridCoupling(debugFsgrid);
   BgBGrid.       setupForGridCoupling(debugFsgrid);
   volGrid.       setupForGridCoupling(debugFsgrid);
   technicalGrid. setupForGridCoupling(debugFsgrid);
   
   // Each dccrg cell may have to communicate with multiple fsgrid cells, if they are on a lower refinement level.
   // Calculate the corresponding fsgrid ids for each dccrg cell and set coupling for each fsgrid id.
   for(auto& dccrgId : cells) {
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);

      for (auto fsgridId : fsgridIds) {

         debugFsgrid = false;
         
         perBGrid.      setGridCoupling(fsgridId, myRank, debugFsgrid);         
         perBDt2Grid.   setGridCoupling(fsgridId, myRank, debugFsgrid);         
         EGrid.         setGridCoupling(fsgridId, myRank, debugFsgrid);
         EDt2Grid.      setGridCoupling(fsgridId, myRank, debugFsgrid);
         EHallGrid.     setGridCoupling(fsgridId, myRank, debugFsgrid);
         EGradPeGrid.   setGridCoupling(fsgridId, myRank, debugFsgrid);
         momentsGrid.   setGridCoupling(fsgridId, myRank, debugFsgrid);
         momentsDt2Grid.setGridCoupling(fsgridId, myRank, debugFsgrid);
         dPerBGrid.     setGridCoupling(fsgridId, myRank, debugFsgrid);
         dMomentsGrid.  setGridCoupling(fsgridId, myRank, debugFsgrid);
         BgBGrid.       setGridCoupling(fsgridId, myRank, debugFsgrid);
         volGrid.       setGridCoupling(fsgridId, myRank, debugFsgrid);
         technicalGrid. setGridCoupling(fsgridId, myRank, debugFsgrid);
      }
   }
   
   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   perBGrid.      finishGridCoupling();
   perBDt2Grid.   finishGridCoupling();
   EGrid.         finishGridCoupling();
   EDt2Grid.      finishGridCoupling();
   EHallGrid.     finishGridCoupling();
   EGradPeGrid.   finishGridCoupling();
   momentsGrid.   finishGridCoupling();
   momentsDt2Grid.finishGridCoupling();
   dPerBGrid.     finishGridCoupling();
   dMomentsGrid.  finishGridCoupling();
   BgBGrid.       finishGridCoupling();
   volGrid.       finishGridCoupling();
   technicalGrid. finishGridCoupling();
   phiprof::stop("Initial fsgrid coupling");

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   // Transfer initial field configuration into the FsGrids
   feedFieldDataIntoFsGrid<fsgrids::N_BFIELD>(mpiGrid,cells,CellParams::PERBX,perBGrid);

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   feedBgFieldsIntoFsGrid(mpiGrid,cells,BgBGrid);
   BgBGrid.updateGhostCells();

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   setupTechnicalFsGrid(mpiGrid, cells, technicalGrid);
   technicalGrid.updateGhostCells();

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   // WARNING this means moments and dt2 moments are the same here.
   feedMomentsIntoFsGrid(mpiGrid, cells, momentsGrid,false);
   feedMomentsIntoFsGrid(mpiGrid, cells, momentsDt2Grid,false);

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   phiprof::start("Init field propagator");
   if (
      initializeFieldPropagator(
         perBGrid,
         perBDt2Grid,
         EGrid,
         EDt2Grid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         volGrid,
         technicalGrid,
         sysBoundaries
      ) == false
   ) {
      logFile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }
   phiprof::stop("Init field propagator");

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__  <<  endl;
   
   // Initialize Poisson solver (if used)
   if (P::propagatePotential == true) {
      phiprof::start("Init Poisson solver");
      if (poisson::initialize(mpiGrid) == false) {
         logFile << "(MAIN): Poisson solver did not initialize correctly!" << endl << writeVerbose;
         exit(1);
      }
      phiprof::stop("Init Poisson solver");
   }

   // Free up memory:
   readparameters.finalize();

   if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   if (P::isRestart == false) {
      // Run Vlasov solver once with zero dt to initialize
      //per-cell dt limits. In restarts, we read the dt from file.
      phiprof::start("compute-dt");
      
      if(P::propagateField) {
         propagateFields(
            perBGrid,
            perBDt2Grid,
            EGrid,
            EDt2Grid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            volGrid,
            technicalGrid,
            sysBoundaries, 0.0, 1.0
         );
      }
      
      calculateSpatialTranslation(mpiGrid,0.0);
      calculateAcceleration(mpiGrid,0.0);
      
      phiprof::stop("compute-dt");
   }
   
   phiprof::start("getVolumeFieldsFromFsGrid");
   // These should be done by initializeFieldPropagator() if the propagation is turned off.
   getVolumeFieldsFromFsGrid(volGrid, mpiGrid, cells);
   phiprof::stop("getVolumeFieldsFromFsGrid");

   if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   // Save restart data
   if (P::writeInitialState) {
      phiprof::start("write-initial-state");
      phiprof::start("fsgrid-coupling-out");
      getFieldDataFromFsGrid<fsgrids::N_BFIELD>(perBGrid,mpiGrid,cells,CellParams::PERBX);
      getFieldDataFromFsGrid<fsgrids::N_EFIELD>(EGrid,mpiGrid,cells,CellParams::EX);
      getFieldDataFromFsGrid<fsgrids::N_EHALL>(EHallGrid,mpiGrid,cells,CellParams::EXHALL_000_100);
      getFieldDataFromFsGrid<fsgrids::N_EGRADPE>(EGradPeGrid,mpiGrid,cells,CellParams::EXGRADPE);
      getDerivativesFromFsGrid(dPerBGrid, dMomentsGrid, BgBGrid, mpiGrid, cells);
      phiprof::stop("fsgrid-coupling-out");
      
      if (myRank == MASTER_RANK)
         logFile << "(IO): Writing initial state to disk, tstep = "  << endl << writeVerbose;
      P::systemWriteDistributionWriteStride.push_back(1);
      P::systemWriteName.push_back("initial-grid");
      P::systemWriteDistributionWriteXlineStride.push_back(0);
      P::systemWriteDistributionWriteYlineStride.push_back(0);
      P::systemWriteDistributionWriteZlineStride.push_back(0);
      P::systemWritePath.push_back("./");

      for(uint si=0; si<P::systemWriteName.size(); si++) {
         P::systemWrites.push_back(0);
      }

      const bool writeGhosts = true;
      if( writeGrid(mpiGrid,&outputReducer,P::systemWriteName.size()-1, writeGhosts) == false ) {
         cerr << "FAILED TO WRITE GRID AT " << __FILE__ << " " << __LINE__ << endl;
      }

      P::systemWriteDistributionWriteStride.pop_back();
      P::systemWriteName.pop_back();
      P::systemWriteDistributionWriteXlineStride.pop_back();
      P::systemWriteDistributionWriteYlineStride.pop_back();
      P::systemWriteDistributionWriteZlineStride.pop_back();
      P::systemWritePath.pop_back();

      phiprof::stop("write-initial-state");
   }

   if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   if (P::isRestart == false) {
      //compute new dt
      phiprof::start("compute-dt");
      getFsGridMaxDt(technicalGrid, mpiGrid, cells);

      if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      computeNewTimeStep(mpiGrid,newDt,dtIsChanged);
      if (P::dynamicTimestep == true && dtIsChanged == true) {
         // Only actually update the timestep if dynamicTimestep is on
         P::dt=newDt;
      }
      phiprof::stop("compute-dt");
   }

   if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;

   if (!P::isRestart) {      
      //go forward by dt/2 in V, initializes leapfrog split. In restarts the
      //the distribution function is already propagated forward in time by dt/2
      phiprof::start("propagate-velocity-space-dt/2");
      if (P::propagateVlasovAcceleration) {
         calculateAcceleration(mpiGrid, 0.5*P::dt);
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }
      phiprof::stop("propagate-velocity-space-dt/2");

   }
   
   phiprof::stop("Initialization");

   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************

   if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
   
   // Main simulation loop:
   if (myRank == MASTER_RANK) logFile << "(MAIN): Starting main simulation loop." << endl << writeVerbose;
   
   phiprof::start("report-memory-consumption");
   report_process_memory_consumption();
   phiprof::stop("report-memory-consumption");
   
   unsigned int computedCells=0;
   unsigned int computedTotalCells=0;
  //Compute here based on time what the file intervals are
   P::systemWrites.clear();
   for(uint i=0;i< P::systemWriteTimeInterval.size();i++){
      int index=(int)(P::t_min/P::systemWriteTimeInterval[i]);
      //if we are already over 1% further than the time interval time that
      //is requested for writing, then jump to next writing index. This is to
      //make sure that at restart we do not write in the middle of
      //the interval.
      if(P::t_min>(index+0.01)*P::systemWriteTimeInterval[i])
         index++;
      P::systemWrites.push_back(index);
   }

   // Invalidate cached cell lists just to be sure (might not be needed)
   P::meshRepartitioned = true;

   unsigned int wallTimeRestartCounter=1;

   int doNow[2]; // 0: writeRestartNow, 1: balanceLoadNow ; declared outside main loop
   int writeRestartNow; // declared outside main loop
   bool overrideRebalanceNow = false; // declared outside main loop
   
   addTimedBarrier("barrier-end-initialization");
   
   phiprof::start("Simulation");
   double startTime=  MPI_Wtime();
   double beforeTime = MPI_Wtime();
   double beforeSimulationTime=P::t_min;
   double beforeStep=P::tstep_min;

   Real nSum = 0.0;
   for(auto cell: cells) {
      creal rho = mpiGrid[cell]->parameters[CellParams::RHOM_R];
      creal dx = mpiGrid[cell]->parameters[CellParams::DX];
      creal dy = mpiGrid[cell]->parameters[CellParams::DY];
      creal dz = mpiGrid[cell]->parameters[CellParams::DZ];
      creal x = mpiGrid[cell]->parameters[CellParams::XCRD];
      creal y = mpiGrid[cell]->parameters[CellParams::YCRD];
      creal z = mpiGrid[cell]->parameters[CellParams::ZCRD];
      
      nSum += rho*dx*dy*dz;
      if(printCells) cout << "Cell " << cell << " rho = " << rho << " x: " << x << " y: " << y << " z: " << z << endl;
   }
   if(printSums) cout << "Rank " << myRank << ", nSum = " << nSum << endl;   
   
   while(P::tstep <= P::tstep_max  &&
         P::t-P::dt <= P::t_max+DT_EPSILON &&
         wallTimeRestartCounter <= P::exitAfterRestarts) {

      if(printLines)  cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      addTimedBarrier("barrier-loop-start");
      
      phiprof::start("IO");

      phiprof::start("checkExternalCommands");
      if(myRank ==  MASTER_RANK) {
         // check whether STOP or KILL or SAVE has been passed, should be done by MASTER_RANK only as it can reset P::bailout_write_restart
         checkExternalCommands();
      }
      phiprof::stop("checkExternalCommands");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      //write out phiprof profiles and logs with a lower interval than normal
      //diagnostic (every 10 diagnostic intervals).
      phiprof::start("logfile-io");
      logFile << "---------- tstep = " << P::tstep << " t = " << P::t <<" dt = " << P::dt << " FS cycles = " << P::fieldSolverSubcycles << " ----------" << endl;
      if (P::diagnosticInterval != 0 &&
          P::tstep % (P::diagnosticInterval*10) == 0 &&
          P::tstep-P::tstep_min >0) {

         phiprof::print(MPI_COMM_WORLD,"phiprof");
         
         double currentTime=MPI_Wtime();
         double timePerStep=double(currentTime  - beforeTime) / (P::tstep-beforeStep);
         double timePerSecond=double(currentTime  - beforeTime) / (P::t-beforeSimulationTime + DT_EPSILON);
         double remainingTime=min(timePerStep*(P::tstep_max-P::tstep),timePerSecond*(P::t_max-P::t));
         time_t finalWallTime=time(NULL)+(time_t)remainingTime; //assume time_t is in seconds, as it is almost always
         struct tm *finalWallTimeInfo=localtime(&finalWallTime);
         logFile << "(TIME) current walltime/step " << timePerStep<< " s" <<endl;
         logFile << "(TIME) current walltime/simusecond " << timePerSecond<<" s" <<endl;
         logFile << "(TIME) Estimated completion time is " <<asctime(finalWallTimeInfo)<<endl;
         //reset before values, we want to report speed since last report of speed.
         beforeTime = MPI_Wtime();
         beforeSimulationTime=P::t;
         beforeStep=P::tstep;
         //report_grid_memory_consumption(mpiGrid);
         report_process_memory_consumption();
      }
      logFile << writeVerbose;
      phiprof::stop("logfile-io");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
// Check whether diagnostic output has to be produced
      if (P::diagnosticInterval != 0 && P::tstep % P::diagnosticInterval == 0) {
         vector<string>::const_iterator it;
         for (it = P::diagnosticVariableList.begin();
              it != P::diagnosticVariableList.end();
         it++) {
            if (*it == "FluxB") {
               phiprof::start("fsgrid-coupling-out");
               getFieldDataFromFsGrid<fsgrids::N_BFIELD>(perBGrid,mpiGrid,cells,CellParams::PERBX);
               phiprof::stop("fsgrid-coupling-out");
            }
            if (*it == "FluxE") {
               phiprof::start("fsgrid-coupling-out");
               getFieldDataFromFsGrid<fsgrids::N_EFIELD>(EGrid,mpiGrid,cells,CellParams::EX);
               phiprof::stop("fsgrid-coupling-out");
            }
         }
         
         phiprof::start("diagnostic-io");
         if (writeDiagnostic(mpiGrid, diagnosticReducer) == false) {
            if(myRank == MASTER_RANK)  cerr << "ERROR with diagnostic computation" << endl;
            
         }
         phiprof::stop("diagnostic-io");
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      bool extractFsGridFields = true;
      // write system, loop through write classes
      for (uint i = 0; i < P::systemWriteTimeInterval.size(); i++) {
         if (P::systemWriteTimeInterval[i] >= 0.0 &&
                 P::t >= P::systemWrites[i] * P::systemWriteTimeInterval[i] - DT_EPSILON) {
            if (extractFsGridFields) {
               vector<string>::const_iterator it;
               for (it = P::outputVariableList.begin();
                    it != P::outputVariableList.end();
               it++) {
                  if (*it == "B" ||
                      *it == "PerturbedB"
                  ) {
                     phiprof::start("fsgrid-coupling-out");
                     getFieldDataFromFsGrid<fsgrids::N_BFIELD>(perBGrid,mpiGrid,cells,CellParams::PERBX);
                     phiprof::stop("fsgrid-coupling-out");
                  }
                  if (*it == "E") {
                     phiprof::start("fsgrid-coupling-out");
                     getFieldDataFromFsGrid<fsgrids::N_EFIELD>(EGrid,mpiGrid,cells,CellParams::EX);
                     phiprof::stop("fsgrid-coupling-out");
                  }
                  if (*it == "HallE") {
                     phiprof::start("fsgrid-coupling-out");
                     getFieldDataFromFsGrid<fsgrids::N_EHALL>(EHallGrid,mpiGrid,cells,CellParams::EXHALL_000_100);
                     phiprof::stop("fsgrid-coupling-out");
                  }
                  if (*it == "GradPeE") {
                     phiprof::start("fsgrid-coupling-out");
                     getFieldDataFromFsGrid<fsgrids::N_EGRADPE>(EGradPeGrid,mpiGrid,cells,CellParams::EXGRADPE);
                     phiprof::stop("fsgrid-coupling-out");
                  }
                  if (*it == "derivs") {
                     phiprof::start("fsgrid-coupling-out");
                     getDerivativesFromFsGrid(dPerBGrid, dMomentsGrid, BgBGrid, mpiGrid, cells);
                     phiprof::stop("fsgrid-coupling-out");
                  }
               }
               extractFsGridFields = false;
            }
            
            phiprof::start("write-system");
            logFile << "(IO): Writing spatial cell and reduced system data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
            const bool writeGhosts = true;
            if( writeGrid(mpiGrid,&outputReducer, i, writeGhosts) == false ) {
               cerr << "FAILED TO WRITE GRID AT" << __FILE__ << " " << __LINE__ << endl;
            }
            P::systemWrites[i]++;
            logFile << "(IO): .... done!" << endl << writeVerbose;
            phiprof::stop("write-system");
         }
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      // Reduce globalflags::bailingOut from all processes
      phiprof::start("Bailout-allreduce");
      MPI_Allreduce(&(globalflags::bailingOut), &(doBailout), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      phiprof::stop("Bailout-allreduce");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      // Write restart data if needed
      // Combined with checking of additional load balancing to have only one collective call.
      phiprof::start("compute-is-restart-written-and-extra-LB");
      if (myRank == MASTER_RANK) {
         if (  (P::saveRestartWalltimeInterval >= 0.0
            && (P::saveRestartWalltimeInterval*wallTimeRestartCounter <=  MPI_Wtime()-initialWtime
               || P::tstep == P::tstep_max
               || P::t >= P::t_max))
            || (doBailout > 0 && P::bailout_write_restart)
            || globalflags::writeRestart
         ) {
            doNow[0] = 1;
            if (globalflags::writeRestart == true) {
               doNow[0] = 2; // Setting to 2 so as to not increment the restart count below.
               globalflags::writeRestart = false; // This flag is only used by MASTER_RANK here and it needs to be reset after a restart write has been issued.
            }
         }
         else {
            doNow[0] = 0;
         }
         if (globalflags::balanceLoad == true) {
            doNow[1] = 1;
            globalflags::balanceLoad = false;
         }
      }
      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      MPI_Bcast( &doNow, 2 , MPI_INT , MASTER_RANK ,MPI_COMM_WORLD);
      writeRestartNow = doNow[0];
      doNow[0] = 0;
      if (doNow[1] == 1) {
         P::prepareForRebalance = true;
         doNow[1] = 0;
      }
      phiprof::stop("compute-is-restart-written-and-extra-LB");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      if (writeRestartNow >= 1){
         phiprof::start("write-restart");
         if (writeRestartNow == 1) {
            wallTimeRestartCounter++;
         }
         
         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         //Write the restart:
         if( writeRestart(mpiGrid,outputReducer,"restart",(uint)P::t, P::restartStripeFactor) == false ) {
            logFile << "(IO): ERROR Failed to write restart!" << endl << writeVerbose;
            cerr << "FAILED TO WRITE RESTART" << endl;
         }
         if (myRank == MASTER_RANK)
            logFile << "(IO): .... done!"<< endl << writeVerbose;
         phiprof::stop("write-restart");
      }
      
      phiprof::stop("IO");
      addTimedBarrier("barrier-end-io");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      //no need to propagate if we are on the final step, we just
      //wanted to make sure all IO is done even for final step
      if(P::tstep == P::tstep_max ||
         P::t >= P::t_max ||
         doBailout > 0) {
         break;
      }
      
      //Re-loadbalance if needed
      //TODO - add LB measure and do LB if it exceeds threshold
      #warning Re-loadbalance has been disabled temporarily for amr debugging
      if(((P::tstep % P::rebalanceInterval == 0 && P::tstep > P::tstep_min) || overrideRebalanceNow) && false) {
         logFile << "(LB): Start load balance, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         balanceLoad(mpiGrid, sysBoundaries);
         addTimedBarrier("barrier-end-load-balance");
         phiprof::start("Shrink_to_fit");
         // * shrink to fit after LB * //
         shrink_to_fit_grid_data(mpiGrid);
         phiprof::stop("Shrink_to_fit");
         logFile << "(LB): ... done!"  << endl << writeVerbose;
         P::prepareForRebalance = false;

         if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
         
         // Re-couple fsgrids to updated grid situation
         phiprof::start("fsgrid-recouple-after-lb");
         
         const vector<CellID>& cells = getLocalCells();

         cout << "Reloadbalance: Local cells are: ";
         for(auto id : cells) cout << id << " ";
         cout << endl;

         debugFsgrid = false;
         
         perBGrid.      setupForGridCoupling(debugFsgrid);
         perBDt2Grid.   setupForGridCoupling(debugFsgrid);        
         EGrid.         setupForGridCoupling(debugFsgrid);
         EDt2Grid.      setupForGridCoupling(debugFsgrid);
         EHallGrid.     setupForGridCoupling(debugFsgrid);
         EGradPeGrid.   setupForGridCoupling(debugFsgrid);
         momentsGrid.   setupForGridCoupling(debugFsgrid);
         momentsDt2Grid.setupForGridCoupling(debugFsgrid);
         dPerBGrid.     setupForGridCoupling(debugFsgrid);
         dMomentsGrid.  setupForGridCoupling(debugFsgrid);
         BgBGrid.       setupForGridCoupling(debugFsgrid);
         volGrid.       setupForGridCoupling(debugFsgrid);
         technicalGrid. setupForGridCoupling(debugFsgrid);

         if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
         
         // Each dccrg cell may have to communicate with multiple fsgrid cells, if they are on a lower refinement level.
         // Calculate the corresponding fsgrid ids for each dccrg cell and set coupling for each fsgrid id.
         // cout << "send tags are: " << endl;
         for(auto& dccrgId : cells) {
            const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
            //cout << "Fsgrid ids for cell " << dccrgId << " are: ";
            for (auto& fsgridId : fsgridIds) {
               //cout << fsgridId << " ";

               debugFsgrid = false;
               
               perBGrid.      setGridCoupling(fsgridId, myRank, debugFsgrid);
               perBDt2Grid.   setGridCoupling(fsgridId, myRank, debugFsgrid);               
               EGrid.         setGridCoupling(fsgridId, myRank, debugFsgrid);
               EDt2Grid.      setGridCoupling(fsgridId, myRank, debugFsgrid);
               EHallGrid.     setGridCoupling(fsgridId, myRank, debugFsgrid);
               EGradPeGrid.   setGridCoupling(fsgridId, myRank, debugFsgrid);
               momentsGrid.   setGridCoupling(fsgridId, myRank, debugFsgrid);
               momentsDt2Grid.setGridCoupling(fsgridId, myRank, debugFsgrid);
               dPerBGrid.     setGridCoupling(fsgridId, myRank, debugFsgrid);
               dMomentsGrid.  setGridCoupling(fsgridId, myRank, debugFsgrid);
               BgBGrid.       setGridCoupling(fsgridId, myRank, debugFsgrid);               
               volGrid.       setGridCoupling(fsgridId, myRank, debugFsgrid);
               technicalGrid. setGridCoupling(fsgridId, myRank, debugFsgrid);
            }
            //cout << endl;
         }
         cout << endl;
         if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
         
         perBGrid.      finishGridCoupling();
         perBDt2Grid.   finishGridCoupling();
         EGrid.         finishGridCoupling();
         EDt2Grid.      finishGridCoupling();
         EHallGrid.     finishGridCoupling();
         EGradPeGrid.   finishGridCoupling();
         momentsGrid.   finishGridCoupling();
         momentsDt2Grid.finishGridCoupling();
         dPerBGrid.     finishGridCoupling();
         dMomentsGrid.  finishGridCoupling();
         BgBGrid.       finishGridCoupling();
         volGrid.       finishGridCoupling();
         technicalGrid. finishGridCoupling();
         phiprof::stop("fsgrid-recouple-after-lb");

         overrideRebalanceNow = false;
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ <<  endl;
      
      //get local cells
      const vector<CellID>& cells = getLocalCells();

      //compute how many spatial cells we solve for this step
      computedCells=0;
      for(size_t i=0; i<cells.size(); i++) {
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            computedCells += mpiGrid[cells[i]]->get_number_of_velocity_blocks(popID)*WID3;
      }
      computedTotalCells+=computedCells;
      
      //Check if dt needs to be changed, and propagate V back a half-step to change dt and set up new situation
      //do not compute new dt on first step (in restarts dt comes from file, otherwise it was initialized before we entered
      //simulation loop
      // FIXME what if dt changes at a restart??
      if(P::dynamicTimestep  && P::tstep > P::tstep_min) {
         getFsGridMaxDt(technicalGrid, mpiGrid, cells);
         computeNewTimeStep(mpiGrid,newDt,dtIsChanged);
         addTimedBarrier("barrier-check-dt");
         if(dtIsChanged) {
            phiprof::start("update-dt");
            //propagate velocity space back to real-time
            if( P::propagateVlasovAcceleration ) {
               // Back half dt to real time, forward by new half dt
               calculateAcceleration(mpiGrid,-0.5*P::dt + 0.5*newDt);
            }
            else {
               //zero step to set up moments _v
               calculateAcceleration(mpiGrid, 0.0);
            }
            
            P::dt=newDt;
            
            logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
            phiprof::stop("update-dt");
            continue; //
            addTimedBarrier("barrier-new-dt-set");
         }
      }
      
      if (P::tstep % P::rebalanceInterval == P::rebalanceInterval-1 || P::prepareForRebalance == true) {
         if(P::prepareForRebalance == true) {
            overrideRebalanceNow = true;
         } else {
            P::prepareForRebalance = true;
         }
         #pragma omp parallel for
         for (size_t c=0; c<cells.size(); ++c) {
            mpiGrid[cells[c]]->get_cell_parameters()[CellParams::LBWEIGHTCOUNTER] = 0;
         }
      }

      phiprof::start("Propagate");
      //Propagate the state of simulation forward in time by dt:
      
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::start("Update system boundaries (Vlasov pre-translation)");
         sysBoundaries.applySysBoundaryVlasovConditions(mpiGrid, P::t+0.5*P::dt); 
         phiprof::stop("Update system boundaries (Vlasov pre-translation)");
         addTimedBarrier("barrier-boundary-conditions");
      }
      
      phiprof::start("Spatial-space");
      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      if( P::propagateVlasovTranslation) {
         calculateSpatialTranslation(mpiGrid,P::dt);
      } else {
         calculateSpatialTranslation(mpiGrid,0.0);
      }

      Real nSum = 0.0;
      for(auto cell: cells) {
         creal rho = mpiGrid[cell]->parameters[CellParams::RHOM_R];
         creal dx = mpiGrid[cell]->parameters[CellParams::DX];
         creal dy = mpiGrid[cell]->parameters[CellParams::DY];
         creal dz = mpiGrid[cell]->parameters[CellParams::DZ];
         creal x = mpiGrid[cell]->parameters[CellParams::XCRD];
         creal y = mpiGrid[cell]->parameters[CellParams::YCRD];
         creal z = mpiGrid[cell]->parameters[CellParams::ZCRD];
         
         nSum += rho*dx*dy*dz;
         if(printCells) cout << "Cell " << cell << " rho = " << rho << " x: " << x << " y: " << y << " z: " << z << endl;
      }
      if(printSums) cout << "Rank " << myRank << ", nSum = " << nSum << endl;   
      
      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      phiprof::stop("Spatial-space",computedCells,"Cells");

      phiprof::start("Compute interp moments");
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHOM_DT2,
         CellParams::VX_DT2,
         CellParams::VY_DT2,
         CellParams::VZ_DT2,
         CellParams::RHOQ_DT2,
         CellParams::P_11_DT2,
         CellParams::P_22_DT2,
         CellParams::P_33_DT2
      );
      phiprof::stop("Compute interp moments");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      // Apply boundary conditions      
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::start("Update system boundaries (Vlasov post-translation)");
         sysBoundaries.applySysBoundaryVlasovConditions(mpiGrid, P::t+0.5*P::dt); 
         phiprof::stop("Update system boundaries (Vlasov post-translation)");
         addTimedBarrier("barrier-boundary-conditions");
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      // Propagate fields forward in time by dt. This needs to be done before the
      // moments for t + dt are computed (field uses t and t+0.5dt)
      if (P::propagateField) {
         phiprof::start("Propagate Fields");

         if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
         
         phiprof::start("fsgrid-coupling-in");
         // Copy moments over into the fsgrid.
         //setupTechnicalFsGrid(mpiGrid, cells, technicalGrid);
         feedMomentsIntoFsGrid(mpiGrid, cells, momentsGrid,false);
         feedMomentsIntoFsGrid(mpiGrid, cells, momentsDt2Grid,true);
         phiprof::stop("fsgrid-coupling-in");

         propagateFields(
            perBGrid,
            perBDt2Grid,
            EGrid,
            EDt2Grid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            volGrid,
            technicalGrid,
            sysBoundaries,
            P::dt,
            P::fieldSolverSubcycles
         );

         phiprof::start("fsgrid-coupling-out");
         // Copy results back from fsgrid.
         getVolumeFieldsFromFsGrid(volGrid, mpiGrid, cells);
         phiprof::stop("fsgrid-coupling-out");
         phiprof::stop("Propagate Fields",cells.size(),"SpatialCells");
         addTimedBarrier("barrier-after-field-solver");
      }

      if (P::propagatePotential == true) {
         poisson::solve(mpiGrid);
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      phiprof::start("Velocity-space");
      if ( P::propagateVlasovAcceleration ) {
         calculateAcceleration(mpiGrid,P::dt);
         addTimedBarrier("barrier-after-ad just-blocks");
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      phiprof::stop("Velocity-space",computedCells,"Cells");
      addTimedBarrier("barrier-after-acceleration");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      phiprof::start("Compute interp moments");
      // *here we compute rho and rho_v for timestep t + dt, so next
      // timestep * //
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHOM,
         CellParams::VX,
         CellParams::VY,
         CellParams::VZ,
         CellParams::RHOQ,
         CellParams::P_11,
         CellParams::P_22,
         CellParams::P_33
      );
      phiprof::stop("Compute interp moments");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
      
      phiprof::stop("Propagate",computedCells,"Cells");
      
      phiprof::start("Project endTimeStep");
      project->hook(hook::END_OF_TIME_STEP, mpiGrid);
      phiprof::stop("Project endTimeStep");

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
    
      // Check timestep
      if (P::dt < P::bailout_min_dt) {
         stringstream s;
         s << "The timestep dt=" << P::dt << " went below bailout.bailout_min_dt (" << to_string(P::bailout_min_dt) << ")." << endl;
         bailout(true, s.str(), __FILE__, __LINE__);
      }
      //Move forward in time
      P::meshRepartitioned = false;
      ++P::tstep;
      P::t += P::dt;

      if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
   }

   double after = MPI_Wtime();

   if(printLines) cout << "I am process " << myRank << " at line " << __LINE__ << " of " << __FILE__ << endl;
   
   phiprof::stop("Simulation");
   phiprof::start("Finalization");
   if (P::propagateField ) { 
      finalizeFieldPropagator();
   }
   if (P::propagatePotential == true) {
      poisson::finalize();
   }
   if (myRank == MASTER_RANK) {
      if (doBailout > 0) {
         logFile << "(BAILOUT): Bailing out, see error log for details." << endl;
      }
      
      double timePerStep;
      if (P::tstep == P::tstep_min) {
         timePerStep=0.0;
      } else {
         timePerStep=double(after  - startTime) / (P::tstep-P::tstep_min);
      }
      double timePerSecond=double(after  - startTime) / (P::t-P::t_min+DT_EPSILON);
      logFile << "(MAIN): All timesteps calculated." << endl;
      logFile << "\t (TIME) total run time " << after - startTime << " s, total simulated time " << P::t -P::t_min<< " s" << endl;
      if(P::t != 0.0) {
         logFile << "\t (TIME) seconds per timestep " << timePerStep  <<
         ", seconds per simulated second " <<  timePerSecond << endl;
      }
      logFile << writeVerbose;
   }
   
   phiprof::stop("Finalization");
   phiprof::stop("main");
   
   phiprof::print(MPI_COMM_WORLD,"phiprof");
   
   if (myRank == MASTER_RANK) logFile << "(MAIN): Exiting." << endl << writeVerbose;
   logFile.close();
   if (P::diagnosticInterval != 0) diagnostic.close();
   
   perBGrid.finalize();
   perBDt2Grid.finalize();
   EGrid.finalize();
   EDt2Grid.finalize();
   EHallGrid.finalize();
   EGradPeGrid.finalize();
   momentsGrid.finalize();
   momentsDt2Grid.finalize();
   dPerBGrid.finalize();
   dMomentsGrid.finalize();
   BgBGrid.finalize();
   volGrid.finalize();
   technicalGrid.finalize();
   
   MPI_Finalize();
   return 0;
}
