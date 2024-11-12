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

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_000_100(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
     (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
     (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
     (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
     (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)-
     (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_010_110(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_001_101(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx ;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_011_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)-
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

// Y
/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_000_010(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_100_110(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)-
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_001_011(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)+
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)+(pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_101_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)-(pC[a_xyy]*pC[a_y])/(12*dy)-
     (pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)+(pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)+
     (pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

// Z
/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_000_001(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)-
     (pC[a_xz]*pC[a_y])/(4*dz)+(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+
     (pC[b_yz]*pC[c_yyz])/(24*dy)-(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)-
     (pC[a_z]*pC[c_xyz])/(24*dx)+(pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_100_101(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       -(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)-
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)-(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_010_011(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       +(pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)
         +(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)-(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_110_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)-(pC[a_xz]*pC[a_y])/(4*dz)-
     (pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)+
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)+
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

// Add helper struct for common calculations
struct HallTermHelpers {
   static Real getHallRhoq(const Real rhoq) {
      return (rhoq <= Parameters::hallMinimumRhoq) ? Parameters::hallMinimumRhoq : rhoq;
   }
   
   static Real getAveragedRhoq(
      const FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
      const int i, const int j, const int k,
      const int i2, const int j2, const int k2) {
      return FOURTH * (
         momentsGrid.get(i, j, k)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i2, j, k)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i, j2, k2)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i2, j2, k2)->at(fsgrids::moments::RHOQ)
      );
   }
};

// Template for calculating Hall term components
template<typename Component>
void calculateEdgeHallTermComponents(
   FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
   FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
   FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
   FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
   FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
   FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
   const FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
   const std::array<Real, Rec::N_REC_COEFFICIENTS>& perturbedCoefficients,
   const int i, const int j, const int k) {
   
   switch(Parameters::ohmHallTerm) {
      case 0:
         return;
         
      case 1: {
         const Real invDx = 1.0 / technicalGrid.DX;
         const Real invDy = 1.0 / technicalGrid.DY;
         const Real invDz = 1.0 / technicalGrid.DZ;
         
         // Direct grid access since we removed const
         const Real Bx = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX) + 
                        BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX);
         const Real By = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBY) + 
                        BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY);
         const Real Bz = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ) + 
                        BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ);
         
         const Real hallRhoq = HallTermHelpers::getHallRhoq(
            momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ));
            
         Component::calculate(EHallGrid, Bx, By, Bz, invDx, invDy, invDz, hallRhoq);
         break;
      }
      
      case 2: {
         // Optimize case 2 with vectorization hints
         #pragma omp simd
         for(int idx = 0; idx < 4; idx++) {
            const Real hallRhoq = HallTermHelpers::getHallRhoq(
               HallTermHelpers::getAveragedRhoq(momentsGrid, i, j, k, i, j, k));
               
            Component::calculateHighOrder(
               EHallGrid, perturbedCoefficients, 
               BgBGrid, technicalGrid, hallRhoq, idx);
         }
         break;
      }
   }
}

// Main Hall term calculation with optimized parallel processing
void calculateHallTermSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool communicateMomentsDerivatives
) {
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0] * gridDims[1] * gridDims[2];

   // Optimize communications
   #pragma omp parallel sections
   {
      #pragma omp section
      dPerBGrid.updateGhostCells();
      
      #pragma omp section
      if(P::ohmGradPeTerm == 0 && communicateMomentsDerivatives) {
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            dMomentsGrid.updateGhostCells();
         } else {
            dMomentsDt2Grid.updateGhostCells();
         }
      }
   }

   // Optimize main computation loops
   #pragma omp parallel for collapse(3) schedule(guided)
   for(FsGridTools::FsIndex_t k = 0; k < gridDims[2]; k++) {
      for(FsGridTools::FsIndex_t j = 0; j < gridDims[1]; j++) {
         for(FsGridTools::FsIndex_t i = 0; i < gridDims[0]; i++) {
            // ... calculate Hall terms using templated functions
            calculateEdgeHallTermComponents<JXBX_000_100>(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
            calculateEdgeHallTermComponents<JXBY_000_010>(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k); 
            calculateEdgeHallTermComponents<JXBZ_000_001>(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
         }
      }
   }
}

// Add helper functions to safely read from const grids
template<typename T, int S>
inline const std::array<T,S>* getGrid(const FsGrid<std::array<T,S>, FS_STENCIL_WIDTH>& grid, int i, int j, int k) {
   // This is safe because we're only reading
   return const_cast<FsGrid<std::array<T,S>, FS_STENCIL_WIDTH>&>(grid).get(i,j,k);
}

template<typename Component>
void calculateEdgeHallTermComponents(
   FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
   FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
   FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
   FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid, 
   FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
   FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
   const FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
   const std::array<Real, Rec::N_REC_COEFFICIENTS>& perturbedCoefficients,
   const int i, const int j, const int k) {
   
   switch(Parameters::ohmHallTerm) {
      case 0:
         return;
         
      case 1: {
         const Real invDx = 1.0 / technicalGrid.DX;
         const Real invDy = 1.0 / technicalGrid.DY;
         const Real invDz = 1.0 / technicalGrid.DZ;
         
         // Use helper function to read from const grids
         const Real Bx = getGrid(perBGrid, i,j,k)->at(fsgrids::bfield::PERBX) + 
                        getGrid(BgBGrid, i,j,k)->at(fsgrids::bgbfield::BGBX);
         const Real By = getGrid(perBGrid, i,j,k)->at(fsgrids::bfield::PERBY) + 
                        getGrid(BgBGrid, i,j,k)->at(fsgrids::bgbfield::BGBY);
         const Real Bz = getGrid(perBGrid, i,j,k)->at(fsgrids::bfield::PERBZ) + 
                        getGrid(BgBGrid, i,j,k)->at(fsgrids::bgbfield::BGBZ);
         
         const Real hallRhoq = HallTermHelpers::getHallRhoq(
            getGrid(momentsGrid, i,j,k)->at(fsgrids::moments::RHOQ));
            
         Component::calculate(EHallGrid, Bx, By, Bz, invDx, invDy, invDz, hallRhoq);
         break;
      }
      
      // ...existing code...
   }
}

// ...rest of the file...
