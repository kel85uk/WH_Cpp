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
#include <../include/boundary_interior_tau_wh.hh>

Matrix<realscalar,4,1> boundary_wh(realscalar x,realscalar Tc,realscalar t,Matrix<realscalar,4,1> qIC,Matrix<realscalar,4,1> wIC,Matrix<realscalar,4,4> S,realscalar L,realscalar l1,realscalar l2,realscalar l3,realscalar l4,realscalar Af,realscalar At,realscalar Cd){
	realscalar eps1 = 1e-15;
	realscalar epse = 1e-3;

	realscalar V0 = qIC(0);
	realscalar P0 = qIC(1);
	realscalar Ux0 = qIC(2);
	realscalar sig0 = qIC(3);

	// Boundary condition matrices
	Matrix<realscalar,4,4> D;
	D << 0,1,0,0,
			 1,0,-1,0,
			 0,0,1,0,
			 0,Af,-Cd,-At;

	// Excitation vector
	Matrix<realscalar,4,1> F;
	F << 0,0,0,0;
	
	Matrix<realscalar,Dynamic,Dynamic> DS;
	Matrix<realscalar,Dynamic,Dynamic> SI;
	DS = D*S;
	SI = S.inverse();

	// Constant coefficients for non-linear valve closure at x = L
	realscalar a1 = SI(0,0);
	realscalar a2 = SI(0,1);
	realscalar a3 = SI(0,2) + a1;
	realscalar a4 = SI(0,3);
	realscalar b1 = SI(2,0);
	realscalar b2 = SI(2,1);
	realscalar b3 = SI(2,2) + b1;
	realscalar b4 = SI(2,3);
	realscalar c1 = D(3,0);
	realscalar c2 = D(3,1);
	realscalar c3 = D(3,2); 
	realscalar c4 = D(3,3);
	realscalar c5 = F(3);
	realscalar d1 = 1;
	realscalar tau = tau_valve(Tc,t);
	realscalar d2 = -std::pow(tau,2)*(V0 - Ux0)*std::abs(V0-Ux0)/P0;
	realscalar d5 = 0;

	// Calculate coefficients alpha, beta
	// For x = 0
	realscalar det13 = DS(0,0)*DS(2,2) - DS(2,0)*DS(0,2);
	realscalar alpha12 =  -(DS(0,1)*DS(2,2) - DS(0,2)*DS(2,1))/det13;
	realscalar alpha14 = -(DS(0,3)*DS(2,2) - DS(0,2)*DS(2,3))/det13;
	realscalar alpha32 = -(DS(2,1)*DS(0,0) - DS(2,0)*DS(0,1))/det13;
	realscalar alpha34 = -(DS(2,3)*DS(0,0) - DS(2,0)*DS(0,3))/det13;

	realscalar beta11 = DS(3,3)/det13;
	realscalar beta31 = -DS(3,1)/det13;
	realscalar beta13 = -DS(1,3)/det13;
	realscalar beta33 = DS(1,1)/det13;
	// For x = L
	realscalar det24 = DS(1,1)*DS(3,3) - DS(3,1)*DS(1,3);
	realscalar alpha21 = -(DS(1,0)*DS(3,3) - DS(1,3)*DS(3,0))/det24;
	realscalar alpha23 = -(DS(1,2)*DS(3,3) - DS(1,3)*DS(3,2))/det24;
	realscalar alpha41 = -(DS(3,0)*DS(1,1) - DS(3,1)*DS(1,0))/det24;
	realscalar alpha43 = -(DS(3,2)*DS(1,1) - DS(3,1)*DS(1,2))/det24;

	realscalar beta22 = DS(3,3)/det24;
	realscalar beta24 = -DS(1,3)/det24;
	realscalar beta42 = -DS(3,1)/det24;
	realscalar beta44 = DS(1,1)/det24;

	realscalar dt3 = L/l3;
	realscalar dt4 = dt3;
	realscalar dt1 = L/l1;
	realscalar dt2 = dt1;
	realscalar zero = 0.;
	Matrix<realscalar,4,1> eta;
	eta = wIC;
//	std::cout.precision(ldbl::digits10);
//	std::cout << "X: " << std::scientific << x-0 << std::endl;
	// Start of boundary computations
	if (t >= eps1){
		if (d_equal_abs_str(x,zero)){
		  eta = boundary_wh(L,Tc,t-dt4,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
		  realscalar eta4 = eta(3);
		  eta = boundary_wh(L,Tc,t-dt2,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
		  realscalar eta2 = eta(1);
		  eta(0) = alpha12*eta2 + alpha14*eta4 + beta11*F(0) + beta13*F(2);
		  eta(1) = eta2;
		  eta(2) = alpha32*eta2 + alpha34*eta4 + beta31*F(0) + beta33*F(2);
		  eta(3) = eta4;
	  }
		else if (d_equal_wk(x,L,epse)){
		  eta = boundary_wh(0,Tc,t-dt3,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
		  realscalar eta3 = eta(2);
		  eta = boundary_wh(0,Tc,t-dt1,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
		  realscalar eta1 = eta(0);
		  // Start nonlinear BC treatment at x = L
		  realscalar a5 = eta1;
		  realscalar b5 = eta3;
		  realscalar det1 = a3*c4-a4*c3;
		  realscalar det2 = b3*c4-b4*c3;
  	    realscalar e1,e2,e5 = 0.;
		  if (d_equal_abs_str(det1,zero) && d_equal_abs_str(det2,zero)){
		    std::cout << "Singular matrix detected for nonlinear BC calculation \n";
	    }
		  else if(d_equal_abs_str(det1,zero) && !d_equal_abs_str(det2,zero)){
//		  (det1 == 0) && (det2 != 0)){
		    e1 = a1*c4 - a4*c1;
		    e2 = a2*c4 - a4*c2;
		    e5 = a5*c4 - a4*c5;
	    }
		  else if(!d_equal_abs_str(det1,zero) && d_equal_abs_str(det2,zero)){
//		  (det1 != 0) && (det2 == 0)){
		    e1 = b1*c4 - b4*c1;
		    e2 = b2*c4 - b4*c2;
		    e5 = b5*c4 - b4*c5;
	    }
		  else if(!d_equal_abs_str(det1,zero) && !d_equal_abs_str(det2,zero)){
		  //(det1 != 0) && (det2 != 0)){
		    e1 = (a1*c4 - a4*c1)/det1 - (b1*c4 - b4*c1)/det2;
		    e2 = (a2*c4 - a4*c2)/det1 - (b2*c4 - b4*c2)/det2;
		    e5 = (a5*c4 - a4*c5)/det1 - (b5*c4 - b4*c5)/det2;
			}
		  realscalar Vr = (d2*e1 + std::sqrt(std::pow(d2*e1,2) - 4*d1*e2*(d2*e5-d5*e2)))/(2*d1*e2);
		  F(1) = F(1) + Vr;// + eta3;
		  
		  eta(1) = alpha21*eta1 + alpha23*eta3 + beta22*F(1) + beta24*F(3);
		  eta(0) = eta1;
		  eta(3) = alpha41*eta1 + alpha43*eta3 + beta42*F(1) + beta44*F(3);
		  eta(2) = eta3;
	  }
		else{
			std::cout << "x is not a boundary point\n";
		}
	}		
	return eta;
}

Matrix<realscalar,4,1> interior_wh(realscalar x,realscalar Tc,realscalar t,Matrix<realscalar,4,1> qIC,Matrix<realscalar,4,1> wIC,Matrix<realscalar,4,4> S,realscalar L,realscalar l1,realscalar l2,realscalar l3,realscalar l4,realscalar Af,realscalar At,realscalar Cd){
	realscalar epe = 1e-15;
	Matrix<realscalar,4,1> eta;
	Matrix<realscalar,4,1> eta3;
	eta = wIC;
	if (t >= epe){
		if ((x>0) && (x<L)){
			eta3 = boundary_wh(0,Tc,t-x/l1,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			realscalar w1 = eta3(0);
			eta3 = boundary_wh(L,Tc,t-(x-L)/l2,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			realscalar w2 = eta3(1);
			eta3 = boundary_wh(0,Tc,t-x/l3,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			realscalar w3 = eta3(2);
			eta3 = boundary_wh(L,Tc,t-(x-L)/l4,qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			realscalar w4 = eta3(3);
			eta << w1,w2,w3,w4;
		}
	}
	return eta;
}

realscalar tau_valve(realscalar Tc,realscalar t){
	realscalar tau = 0;
	if(t >= 0 && t <= 0.4*Tc)
		tau = std::pow(1-t/Tc,(realscalar)3.53);
	else if(t > 0.4*Tc && t <= Tc)
		tau = 0.394*std::pow(1-t/Tc,(realscalar)1.70);
	return tau;
}
