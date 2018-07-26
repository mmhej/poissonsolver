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
	const double r0 = 1.0;
	const double m = 4.0; 

	const double domain_xmin[2]   = {-2.0, -2.0};
	const double domain_xmax[2]   = { 2.0,  2.0};

	int domain_bounds[2] = {0,0};

//	int domain_ncell[2]  = { 64, 64};
	int domain_ncell[2]  = { 128, 128};
//	int domain_ncell[2]  = { 256, 256};
//	int domain_ncell[2]  = { 512, 512};
//	int domain_ncell[2]  = { 1024, 1024};
//	int domain_ncell[2]  = { 2048, 2048};
//	int domain_ncell[2]  = { 3000, 3000};


//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     d;
	int     i, j, jn, ij, ji, pq;
	int     n;
	double  x, y, r;
	int     ncell[2];
	double  dx[2], xmin[2], xmax[2];

	double err, error;
	double nrm, norm;
	double solX, solY, solZ;

	double * Ax;
	double * Ay;
	double * Az;

	double * Bx;
	double * By;
	double * Bz;

	double diffX, diffY;

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
	green.setup2d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	ncell[0] = green.partition[rank].ncell[0];
	ncell[1] = green.partition[rank].ncell[1];

	xmin[0]  = domain_xmin[0] + dx[0] * double(green.partition[rank].icell[0]);
	xmin[1]  = domain_xmin[1] + dx[1] * double(green.partition[rank].icell[1]);

//----------------------------------------------------------------------------//
// Allocate fields
//----------------------------------------------------------------------------//
	Ax = new double[ncell[0] * ncell[1]]();
	Ay = new double[ncell[0] * ncell[1]]();
//	Az = new double[ncell[0] * ncell[1]]();
	Bx = new double[ncell[0] * ncell[1]]();
//	By = new double[ncell[0] * ncell[1]]();
//	Bz = new double[ncell[0] * ncell[1]]();

//----------------------------------------------------------------------------//
// Initiate fields
//----------------------------------------------------------------------------//
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
/*
				Bx[ij] = 4.0 * c * pow(r0,2) 
				       * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
				       * ( pow(r0,4)
				         - pow(x,4)
				         - pow(y,4)
				         - 2.0*pow(x,2)*pow(y,2)
				         - c*pow(x,2)*pow(r0,2)
				         - c*pow(y,2)*pow(r0,2) )
				       * pow(pow(r0,2) - pow(x,2) - pow(y,2),-4);
*/

				Bx[ij] = pow(1.0 - r*r, m);

			}
			else
			{
				Bx[ij] = 0.0;
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

	green.push( Bx, NULL, NULL, NULL, NULL, NULL);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Solve
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Solve" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	green.solve2d( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Map back to ClientArray
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [greenfish]: Pull" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

//	green.pull( B, NULL, NULL, A, NULL, NULL );
	green.pull( Bx, NULL, NULL, Ax, Ay, NULL );

//----------------------------------------------------------------------------//
// Calculate error integral
//----------------------------------------------------------------------------//
	err   = 0.0;
	error = 0.0;
	nrm   = 0.0;
	norm  = 0.0;

	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
		y  = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x  = xmin[0] + (double(i) + 0.5)*dx[0];

			r = sqrt(x*x + y*y);


			if(r < 0.25*dx[0])
			{

			}
			else if( r < r0 )
			{
/*
				solX = - 2.0 * c * pow(r0,2) * x
				     * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
				     * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);

				solY = - 2.0 * c * pow(r0,2) * y
				     * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
				     * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);
*/
				solX = - x * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2));
				solY = - y * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2));

				err += dx[0]*dx[1]* pow( Ax[ij] - solX , 2);
				nrm += dx[0]*dx[1]* pow( solX , 2);

				err += dx[0]*dx[1]* pow( Ay[ij] - solY , 2);
				nrm += dx[0]*dx[1]* pow( solY , 2);
			}
			else
			{

				solX = - x/(2.0*(m+1.0)*pow(r,2));
				solY = - y/(2.0*(m+1.0)*pow(r,2));

				err += dx[0]*dx[1]* pow( Ax[ij] - solX , 2);
				nrm += dx[0]*dx[1]* pow( solX , 2);

				err += dx[0]*dx[1]* pow( Ay[ij] - solY , 2);
				nrm += dx[0]*dx[1]* pow( solY , 2);


/*
				err += dx[0]*dx[1]* pow( Ax[ij] , 2);
				err += dx[0]*dx[1]* pow( Ay[ij] , 2);
*/
			}


		}
	}

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
	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
		y  = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x  = xmin[0] + (double(i) + 0.5)*dx[0];

			r = sqrt(x*x + y*y);


			if(r < 0.25*dx[0])
			{
				solX = 0.0;
				solY = 0.0;

				diffX = Ax[ij] - solX;
				diffY = Ay[ij] - solY;

//				diffX = solX;
//				diffY = solY;

			}
			else if( r < r0 )
			{
// A
//				err = Ax[ij] - exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2)));

// B
//				err = Bx[ij] - (4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4));


				solX = - x * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2));
				solY = - y * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2));

				diffX = Ax[ij] - solX;
				diffY = Ay[ij] - solY;


//				diffX = solX;
//				diffY = solY;

// dAdX
//				diffX = Ax[ij] - (- 2.0 * c * pow(r0,2) * x * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) );

// dAdY
//				diffY = Ay[ij] - (- 2.0 * c * pow(r0,2) * y * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2) );

			}
			else
			{
//				err = Bx[ij];

				solX = - x/(2.0*(m+1.0)*pow(r,2));
				solY = - y/(2.0*(m+1.0)*pow(r,2));

				diffX = Ax[ij] - solX;
				diffY = Ay[ij] - solY;

//				diffX = solX;
//				diffY = solY;

			}

			err = sqrt( pow(diffX,2) + pow(diffY,2) );

			outfile << std::scientific << std::setw(17) << x << std::setw(17) << y << std::setw(17) << err << "\n";

		}
	}
	outfile.close();



//----------------------------------------------------------------------------//
// Finalize OpenMPI
//----------------------------------------------------------------------------//
	MPI_Finalize();

	return 0;
}

