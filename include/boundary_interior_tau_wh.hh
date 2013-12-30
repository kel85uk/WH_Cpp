/*    Copyright (C) 2013  kklloh

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#ifndef boundary_interior_tau_HH
#define boundary_interior_tau_HH
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <../include/Types.hh>
#include <../include/double_compares.h>

using namespace Eigen;
using namespace double_compares;

Matrix<realscalar,4,1> boundary_wh(realscalar x,realscalar Tc,realscalar t,Matrix<realscalar,4,1> qIC,Matrix<realscalar,4,1> wIC,Matrix<realscalar,4,4> S,realscalar L,realscalar l1,realscalar l2,realscalar l3,realscalar l4,realscalar Af,realscalar At,realscalar Cd);

Matrix<realscalar,4,1> interior_wh(realscalar x,realscalar Tc,realscalar t,Matrix<realscalar,4,1> qIC,Matrix<realscalar,4,1> wIC,Matrix<realscalar,4,4> S,realscalar L,realscalar l1,realscalar l2,realscalar l3,realscalar l4,realscalar Af,realscalar At,realscalar Cd);

realscalar tau_valve(realscalar Tc,realscalar t);

#endif
