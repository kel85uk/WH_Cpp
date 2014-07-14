/** Full FSI version - Main Program for Delft Hydraulics Problem
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 kklloh
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include "../include/Eigen/Dense"
#include "../include/Eigen/LU"
#include "../include/boundary_interior_tau_wh.hh"
#include "../include/Types.hh"
#include "mpi.h"

#ifndef PI
#define PI (4*atan(1))
#endif

#define WORKTAG 1
#define DIETAG 2

using namespace Eigen;

void Masterprocess(const realscalar&);
void Slaveprocess();
std::string convertInt(int number);

int numprocs, myrank;
MPI_Status status;

std::string convertReal(const realscalar &number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

int main(int argc, char *argv[]){
	int result = MPI_Init(&argc,&argv);
	double start, end;
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Barrier(MPI_COMM_WORLD);
	realscalar cdd = std::atof(argv[1]);
	start = MPI_Wtime();
	if(myrank == 0){
		Masterprocess(cdd);
//		std::cout << "Ok! \n" << std::endl;
	}
	else{
		Slaveprocess();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	MPI_Finalize();
	if (myrank == 0){ /* use time on master node */
	std::string timings_name;
	timings_name = std::string("timings_");
	timings_name += convertInt(numprocs);
	timings_name += std::string(".dat");
	std::ofstream file_timings(timings_name);
	if (file_timings.is_open()){
	  file_timings << end-start << '\n';
	}
//    std::cout << "Runtime = " << end-start << " sec" << std::endl;
  }
	return 0;
}

void Masterprocess(const realscalar& cdd){
		realscalar L = 20.;
		realscalar R = 0.3985;//0.3587; //398.5e-3;
		realscalar ee = 0.008;//0.001;//8e-3;
		realscalar E = 210e9;
		realscalar rho_t = 7985;
		realscalar nu = 0.29;
		realscalar K_w = 2.14e9;
		realscalar rho_f = 998;
		realscalar Tc = 0.03;//0.003;//0.03;
		realscalar Q0 = 0.5;
		realscalar zeta0 = 0.2;
		realscalar Cd = cdd;//5.6104e18;

		realscalar Af = PI*std::pow(R,2);
		realscalar At = PI*(std::pow(R+ee,2) - std::pow(R,2));

		realscalar c2_f = 1/(rho_f*(1/K_w + 2*R*(1-std::pow(nu,2))/(E*ee)));
		realscalar c2_t = E/(rho_t);

		realscalar c_f = std::sqrt(c2_f);
		realscalar c_t = std::sqrt(c2_t);
		realscalar g_2 = c2_f + c2_t + 2*std::pow(nu,2)*(rho_f/rho_t)*R*c2_f/ee;

		realscalar aa = 1;
		realscalar bb = -g_2;
		realscalar cc = c2_f*c2_t;

		realscalar l2_1 = (-bb-sqrt(std::pow(bb,2) - 4*aa*cc))/(2*aa);
		realscalar l2_3 = (-bb+sqrt(std::pow(bb,2) - 4*aa*cc))/(2*aa);
		realscalar l2_2 = l2_1;
		realscalar l2_4 = l2_3;

		realscalar l1 = std::sqrt(l2_1);
		realscalar l2 = -l1;
		realscalar l3 = std::sqrt(l2_3);
		realscalar l4 = -l3;

		Matrix<realscalar,4,4> A;
		A << 1,0,0,0,
				 0,1/(rho_f*c2_f),0,0,
				 0,0,1,0,
				 0,nu*R/(E*ee),0,-1/(rho_t*c2_t);
		Matrix<realscalar,4,4> B;
		B << 0,1/rho_f,0,0,
				 1,0,0,0,
				 0,0,0,-1/rho_t,
				 0,0,1,0;

		realscalar al1 = c2_t/l2_1 - 1;
		realscalar al2 = c2_t/l2_2 - 1;
		realscalar t11 = 1;
		realscalar t21 = 1;
		realscalar t12 = l1;
		realscalar t22 = l2;
		realscalar t13 = 2*nu/al1;
		realscalar t23 = 2*nu/al2;
		realscalar t14 = c2_t*t13/l1;
		realscalar t24 = c2_t*t23/l2;
		realscalar al3 = c2_f/l2_3 - 1;
		realscalar al4 = c2_f/l2_4 - 1;
		realscalar t34 = l3;
		realscalar t44 = l4;
		realscalar t33 = l2_3/c2_t;
		realscalar t43 = l2_4/c2_t;
		realscalar t32 = rho_f*nu/E*R/ee*c2_f*l3/al3;
		realscalar t42 = rho_f*nu/E*R/ee*c2_f*l4/al4;
		realscalar t31a = rho_f*nu/E*R/ee*l2_3 + l3*t32/c2_f;
		realscalar t41a = rho_f*nu/E*R/ee*l2_4 + l4*t42/c2_f;
		realscalar t31 = t32/l3;
		realscalar t41 = t42/l4;
		Matrix<realscalar,4,4> T;
		T << t11,t12,t13,t14,
				 t21,t22,t23,t24,
				 t31,t32,t33,t34,
				 t41,t42,t43,t44;
		Matrix<realscalar,Dynamic,Dynamic> S;
		S = (T*A).inverse();

		//Define initial conditions
		realscalar V0 = Q0/Af;
		realscalar P0 = zeta0*1/2*rho_f*V0*std::abs(V0);
		realscalar Pj0 = rho_f*V0*l1;
		Matrix<realscalar,4,1> qIC;
		qIC	<< V0,P0,0,Af*P0/At;
		Matrix<realscalar,Dynamic,1> wIC;
		wIC = S.inverse()*qIC;

		// Test case
		realscalar xtest = 10;
		realscalar Tend = 0.2;
		realscalar dt = 0.001;// %1.45673544858331e-4;
		unsigned int N = std::ceil(Tend/dt);
		Matrix<realscalar,Dynamic,1> ti(N+1);
		for (int i = 0; i < N+1; ++i)
			ti(i) = (i)*dt;
		int received_answers = 0;
		int who,whomax;
		Matrix<realscalar,Dynamic,Dynamic> W0(4,N+1);
		Matrix<realscalar,Dynamic,Dynamic> WL(4,N+1);
		Matrix<realscalar,Dynamic,Dynamic> q0vec(4,N+1);
		Matrix<realscalar,Dynamic,Dynamic> qLvec(4,N+1);	
		Matrix<realscalar,Dynamic,Dynamic> wtest(4,N+1);
		Matrix<realscalar,Dynamic,Dynamic> qtestvec(4,N+1);
		Matrix<realscalar,Dynamic,1> Sendrecvpack(13);
		// Broadcast the initial conditions and constants to the other processors
		MPI_Bcast(qIC.data(),qIC.size(),MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(wIC.data(),wIC.size(),MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&N,1,MPI_INT,myrank,MPI_COMM_WORLD);
		MPI_Bcast(ti.data(),ti.size(),MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(S.data(),S.size(),MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&L,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&l1,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&l2,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&l3,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&l4,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&Af,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&At,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&Cd,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&Tc,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);
		MPI_Bcast(&xtest,1,MPI_DOUBLE,myrank,MPI_COMM_WORLD);

		int total_work = N+1;
		MPI_Request request;
		int i, ii = 0; // We haven't finished any useful work
		// Send task to slaves
		whomax = numprocs - 1;
		if (whomax > total_work)
			whomax = total_work;
		for (who = 1; who <= whomax; ++who, ++ii){
			MPI_Send(&ii,1,MPI_INT,who,WORKTAG,MPI_COMM_WORLD);
		}
//		for (int i = 0; i < N+1; ++i){
		while( received_answers < total_work ){
			MPI_Recv(Sendrecvpack.data(),Sendrecvpack.size(),MPI_DOUBLE,MPI_ANY_SOURCE,WORKTAG,MPI_COMM_WORLD,&status);
			who = status.MPI_SOURCE;
			i = std::round(Sendrecvpack(12));
//			std::cout << "Receive datapack " << i << "from process " << who << std::endl;
			W0.col(i) = Sendrecvpack.segment(0,4);//boundary_wh(0.,Tc,ti(i),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			WL.col(i) = Sendrecvpack.segment(4,4);//boundary_wh(L,Tc,ti(i),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			q0vec.col(i) = S*W0.col(i);
			qLvec.col(i) = S*WL.col(i);
			wtest.col(i) = Sendrecvpack.segment(8,4);//interior_wh(xtest,Tc,ti(i),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			qtestvec.col(i) = S*wtest.col(i);
//			std::cout << Sendrecvpack(12) << std::endl;
			++received_answers;
			if (ii < N+1){
				MPI_Send(&ii,1,MPI_INT,who,WORKTAG,MPI_COMM_WORLD);
//				std::cout << "Waiting to receive results" << std::endl;
				++ii;
			}
		}
		for (who = 1; who < numprocs; ++who){
			MPI_Send(&xtest,1,MPI_INT,who,DIETAG,MPI_COMM_WORLD);
		}
		
		std::string timename, qLname, q0name, qtestname;
		std::string numberprocs(convertInt(numprocs));
		std::string scd(convertReal(Cd));
		timename = std::string("time_") + numberprocs + std::string("_") + scd + std::string(".dat");
		qLname = std::string("qL_") + numberprocs + std::string("_") + scd + std::string(".dat");
		q0name = std::string("q0_") + numberprocs + std::string("_") + scd + std::string(".dat");
		qtestname = std::string("qtest_") + numberprocs + std::string("_") + scd + std::string(".dat");
		std::ofstream file_t(timename.c_str());
		if (file_t.is_open()){
		  file_t << ti << '\n';
		}
		std::ofstream file_qL(qLname.c_str());
		if (file_qL.is_open()){
		  file_qL << qLvec.transpose() << '\n';
		}
		std::ofstream file_q0(q0name.c_str());
		if (file_q0.is_open()){
		  file_q0 << q0vec.transpose() << '\n';
		}
		std::ofstream file_qtest(qtestname.c_str());
		if (file_qtest.is_open()){
		  file_qtest << qtestvec.transpose() << '\n';
		}
		//Matrix<realscalar,Dynamic,1> Pm(3);
		//Pm	<< qLvec.cwiseMax(2),q0vec.cwiseMax(2),qtestvec.cwiseMax(2);
		Matrix<realscalar,Dynamic,Dynamic> qLT, q0T, qtestT;
		qLT = qLvec.transpose();
		q0T = q0vec.transpose();
		qtestT = qtestvec.transpose();
		Matrix<realscalar,3,1> PMvec;
//		PMvec << qLT.col(1).maxCoeff(),q0T.col(1).maxCoeff(),qtestT.col(1).maxCoeff();
		realscalar sss1 = qLT.col(1).maxCoeff(), sss2 = q0T.col(1).maxCoeff(), sss3 = qtestT.col(1).maxCoeff();
		PMvec << sss1,sss2,sss3;
		std::cout << cdd << "\t" << std::setprecision(12) << PMvec.maxCoeff()/Pj0 << "\t" << Pj0 << std::endl;
//		std::cout << "max|P|/Pj = " << ratP << std::endl;
}

void Slaveprocess(){
	Matrix<realscalar,4,4> S;
	Matrix<realscalar,4,1> qIC;
	Matrix<realscalar,4,1> wIC;
	realscalar L,l1,l2,l3,l4,Af,At,Cd,Tc,xtest;
	int tag,N,job;
	MPI_Bcast(qIC.data(),qIC.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(wIC.data(),wIC.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
	Matrix<realscalar,Dynamic,1> ti(N+1);
	MPI_Bcast(ti.data(),ti.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(S.data(),S.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&l1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&l2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&l3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&l4,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Af,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&At,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Cd,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Tc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&xtest,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	Matrix<realscalar,Dynamic,1> Sendrecvpack(13);
	Matrix<realscalar,4,1> W0;
	Matrix<realscalar,4,1> WL;
	Matrix<realscalar,4,1> Wtest;
	MPI_Recv(&job,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	if (status.MPI_TAG != DIETAG){
		while(true){
			W0 = boundary_wh(0.,Tc,ti(job),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			WL = boundary_wh(L,Tc,ti(job),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
			Wtest = interior_wh(xtest,Tc,ti(job),qIC,wIC,S,L,l1,l2,l3,l4,Af,At,Cd);
//			std::cout << "Packing data from processor " << myrank << std::endl;
			Sendrecvpack.segment(0,4) = W0;
//			std::cout << "Packing data 1 from processor " << myrank << std::endl;
			Sendrecvpack.segment(4,4) = WL;
//			std::cout << "Packing data 2 from processor " << myrank << std::endl;
			Sendrecvpack.segment(8,4) = Wtest;
//			std::cout << "Packing data 3 from processor " << myrank << std::endl;
			Sendrecvpack(12) = (double)job;
			MPI_Send(Sendrecvpack.data(),Sendrecvpack.size(),MPI_DOUBLE,0,WORKTAG,MPI_COMM_WORLD);
//			std::cout << "Sent data to root\n";
			MPI_Recv(&job,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);	
			if(status.MPI_TAG == DIETAG)
				break;
		}
	}
}

std::string convertInt(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
