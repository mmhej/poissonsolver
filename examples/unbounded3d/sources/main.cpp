//----------------------------------------------------------------------------//
/*
  File:         main.cpp

  Description:  
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// INCLUDE HEADERS
//----------------------------------------------------------------------------//
// System
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
// External
#include "mpi.h"
#include "class_greenfish.hpp"
// Internal

//----------------------------------------------------------------------------//
// MAIN PROGRAM
//----------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

	const double c = 10.0;
	const double r0 = 0.5;

	const double domain_xmin[3]   = { -0.5, -1.0, -1.0 };
	const double domain_xmax[3]   = {  0.5,  1.0,  1.0 };

	int domain_bounds[3] = { 0, 0, 0 };

	int domain_ncell[3]  = { 64, 64, 64 };

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     d;
	int     i, j, k, kn, kjn, ijk;
	int     n;
	double  x, y, z, rho, phi, theta;
	int     ncell[3];
	double  dx[3], xmin[3], xmax[3];

	double err, error;
	double nrm, norm;
	double Bmag;

	double * Ax;
	double * Ay;
	double * Az;

	double * Bx;
	double * By;
	double * Bz;

	double * Vx;
	double * Vy;
	double * Vz;

	double diffX, diffY, diffZ;

	std::ostringstream str;
	std::string filename;

//----------------------------------------------------------------------------//
// Initialize the OpenMPI library
//----------------------------------------------------------------------------//
	MPI_Init(&argc, &argv);

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{
		std::cout << std::endl;
		std::cout << "INITIATING TEST... " << std::endl;
		std::cout << std::endl;
		std::cout << "Number of processors intiated: " << nproc << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

//	for( d = 0; d < 5; ++d) // convergence test loop
//	{
//	domain_ncell[0] = pow(2, d+4);
//	domain_ncell[1] = pow(2, d+4);

//----------------------------------------------------------------------------//
// Setup domain
//----------------------------------------------------------------------------//
	dx[0] = ( domain_xmax[0] - domain_xmin[0] )/double(domain_ncell[0]);
	dx[1] = ( domain_xmax[1] - domain_xmin[1] )/double(domain_ncell[1]);
	dx[2] = ( domain_xmax[2] - domain_xmin[2] )/double(domain_ncell[2]);

//----------------------------------------------------------------------------//
// Setup GreenFish
//----------------------------------------------------------------------------//
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Setup" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	class_greenfish green;
	green.lhs_grad = true; // specify lhs operator
	green.setup3d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	ncell[0] = green.partition[rank].ncell[0];
	ncell[1] = green.partition[rank].ncell[1];
	ncell[2] = green.partition[rank].ncell[2];

	xmin[0]  = domain_xmin[0] + dx[0] * double(green.partition[rank].icell[0]);
	xmin[1]  = domain_xmin[1] + dx[1] * double(green.partition[rank].icell[1]);
	xmin[2]  = domain_xmin[2] + dx[2] * double(green.partition[rank].icell[2]);

//----------------------------------------------------------------------------//
// Allocate fields
//----------------------------------------------------------------------------//
	Ax    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Ay    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Az    = new double[ncell[0] * ncell[1] * ncell[2]]();

	Bx    = new double[ncell[0] * ncell[1] * ncell[2]]();
	By    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Bz    = new double[ncell[0] * ncell[1] * ncell[2]]();

//	dAdX = new double[ncell[0] * ncell[1]]();
//	dAdY = new double[ncell[0] * ncell[1]]();

//----------------------------------------------------------------------------//
// Initiate fields
//----------------------------------------------------------------------------//
	for (k = 0; k < ncell[2]; ++k )
	{
		kn = k * ncell[1];
		z  = xmin[2] + (double(k) + 0.5)*dx[2];
		for (j = 0; j < ncell[1]; ++j )
		{
			kjn = (kn + j) * ncell[0];
			y  = xmin[1] + (double(j) + 0.5)*dx[1];
			for (i = 0; i < ncell[0]; ++i )
			{
				ijk = kjn + i;
				x  = xmin[0] + (double(i) + 0.5)*dx[0];

				rho   = sqrt(z*z + y*y);
				phi   = sqrt( pow(rho - r0,2) + pow(x,2) );
				theta = atan2(z,y);

				if( phi < r0)
				{
					Bmag = -exp(- c*pow(r0,2)/(2.0*r0*rho - pow(rho,2) - pow(x,2))) *
					     ( 4.0*pow(c,2)*pow(r0,4)*pow(x,2)*pow(rho,2)
					     - 16.0*pow(r0,4)*pow(rho,4)
					     + 32.0*pow(r0,3)*pow(rho,5)
					     - 24.0*pow(r0,2)*pow(rho,6)
					     + 8.0*r0*pow(rho,7) 
					     - 4.0*pow(rho,6)*pow(x,2)
					     - 6.0*pow(rho,4)*pow(x,4) 
					     - 4.0*pow(rho,2)*pow(x,6)
					     - 8.0*c*pow(r0,5)*pow(rho,3)
					     + 8.0*c*pow(r0,4)*pow(rho,4)
					     - 6.0*c*pow(r0,3)*pow(rho,5) 
					     + 4.0*pow(c,2)*pow(r0,6)*pow(rho,2)
					     - 8.0*pow(c,2)*pow(r0,5)*pow(rho,3) 
					     + 4.0*pow(c,2)*pow(r0,4)*pow(rho,4) 
					     + 2.0*c*pow(r0,2)*pow(rho,6)
					     + 32.0*pow(r0,3)*pow(rho,3)*pow(x,2)
					     - 48.0*pow(r0,2)*pow(rho,4)*pow(x,2)
					     - 24.0*pow(r0,2)*pow(rho,2)*pow(x,4) 
					     + 24.0*r0*pow(rho,5)*pow(x,2)
					     + 24.0*r0*pow(rho,3)*pow(x,4)
					     + 8.0*r0*rho*pow(x,6)
					     + 2.0*c*pow(r0,3)*rho*pow(x,4)
					     + 2.0*c*pow(r0,2)*pow(rho,2)*pow(x,4)
					     - 4.0*c*pow(r0,3)*pow(rho,3)*pow(x,2) 
					     + 4.0*c*pow(r0,2)*pow(rho,4)*pow(x,2) 
					     - pow(rho,8) - pow(x,8))
					     * pow(2.0*r0*rho - pow(rho,2) - pow(x,2),-4) * pow(rho,-2);

					Bx[ijk] =   0.0;
					By[ijk] = - sin(theta)*Bmag;
					Bz[ijk] =   cos(theta)*Bmag;
				}
				else
				{
					Bx[ijk] = 0.0;
					By[ijk] = 0.0;
					Bz[ijk] = 0.0;
				}

			}
		}
	}

//----------------------------------------------------------------------------//
// Solve Poisson equation
//----------------------------------------------------------------------------//
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Push to GreenFish array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Push" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

//	green.push( B, NULL, NULL, NULL, NULL, NULL);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Solve
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Solve" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

//	green.solve2d( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Map back to ClientArray
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Pull" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

//	green.pull( B, NULL, NULL, A, NULL, NULL );
//	green.pull( B, NULL, NULL, dAdX, dAdY, NULL );

//----------------------------------------------------------------------------//
// Calculate error integral
//----------------------------------------------------------------------------//
	err   = 0.0;
	error = 0.0;
	nrm   = 0.0;
	norm  = 0.0;
/*
	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
		y  = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x  = xmin[0] + (double(i) + 0.5)*dx[0];

			r = sqrt(x*x + y*y);
			if( r < r0 )
			{
// A
//				err += dx[0]*dx[1]* pow( A[ij] - exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2))) , 2);
//				nrm += dx[0]*dx[1]* pow( exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2))), 2 );

// B
//				err += dx[0]*dx[1]* pow( bump[ij] - (4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4)) ,2);
//				nrm += dx[0]*dx[1]* pow( 4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4) ,2);

// dAdX
				err += dx[0]*dx[1]* pow( dAdX[ij] - (- 2.0 * c * pow(r0,2) * x * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) ) , 2);
				nrm += dx[0]*dx[1]* pow( - 2.0 * c * pow(r0,2) * y * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) , 2);

// dAdY
				err += dx[0]*dx[1]* pow( dAdY[ij] - (- 2.0 * c * pow(r0,2) * y * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) ) , 2);
				nrm += dx[0]*dx[1]* pow( - 2.0 * c * pow(r0,2) * x * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) , 2);

			}
			else
			{

//				err += dx[0]*dx[1]* pow( A[ij], 2);
//				err += dx[0]*dx[1]* pow( B[ij], 2);
				err += dx[0]*dx[1]* pow( dAdX[ij] , 2);
				err += dx[0]*dx[1]* pow( dAdY[ij] , 2);
			}


		}
	}
*/
	MPI_Reduce( &err, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce( &nrm, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	if(rank == 0)
	{
		std::cout << "Error: " << std::scientific << std::setw(7) << ncell[0] << std::setw(17) << dx[0] << std::setw(17) << sqrt( error/norm ) << std::endl;
	}

//	} // convergence test loop


//----------------------------------------------------------------------------//
// Output field
//----------------------------------------------------------------------------//
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << "  Output fields" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif
/*
	str << std::setw(2) << std::setfill('0') << rank;
	filename = "./output/mesh_P" + str.str();
	str.str(""); // clear str
	std::ofstream outfile( filename.c_str() );
//	outfile.precision(16);
	outfile.precision(8);
	if(!outfile.is_open())
	{
		std::cerr << "ERROR: cannot open outfile." << std::endl;
		return 0;
	}



	outfile.close();
*/


//----------------------------------------------------------------------------//
// Finalize OpenMPI
//----------------------------------------------------------------------------//
	MPI_Finalize();

	return 0;
}

