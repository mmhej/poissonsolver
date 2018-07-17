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

	const double domain_xmin[2]   = { 0.0,  0.0};
	const double domain_xmax[2]   = { 1.0,  1.0};

	int domain_bounds[2] = {1,1};

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
	int     i, j, ij, ji, pq;
	int     n;
	double  x, y, r;
	int     ncell[2];
	double  dx[2], xmin[2], xmax[2];

	double err, error;
	double nrm, norm;

	double * A;
	double * B;
	double * vX;
	double * vY;

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
//	green.lhs_grad = true; // specify lhs operator
	green.setup2d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	ncell[0] = green.partition.real[rank].ncell[0];
	ncell[1] = green.partition.real[rank].ncell[1];

	xmin[0]  = domain_xmin[0] + dx[0] * double(green.partition.real[rank].icell[0]);
	xmin[1]  = domain_xmin[1] + dx[1] * double(green.partition.real[rank].icell[1]);

//----------------------------------------------------------------------------//
// Allocate fields
//----------------------------------------------------------------------------//
	A  = new double[ncell[0] * ncell[1]]();
	B  = new double[ncell[0] * ncell[1]]();
	vX = new double[ncell[0] * ncell[1]]();
	vY = new double[ncell[0] * ncell[1]]();

//----------------------------------------------------------------------------//
// Initiate fields
//----------------------------------------------------------------------------//
	for (j = 0; j < ncell[1]; ++j )
	{
		y = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			x = xmin[0] + (double(i) + 0.5)*dx[0];
			ij = j * ncell[0] + i;

			B[ij] = 8.0 * pi * pi * sin(2.0*pi*x) * sin(2.0*pi*y);

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

	green.push2d( NULL, NULL, B, NULL, NULL, NULL);

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

	A  = new double[ncell[0] * ncell[1]]();
	vX = new double[ncell[0] * ncell[1]]();
	vY = new double[ncell[0] * ncell[1]]();
	green.pull2d( NULL, NULL, B, vX, vY, A );

//----------------------------------------------------------------------------//
// Calculate error
//----------------------------------------------------------------------------//
	err   = 0.0;
	error = 0.0;
	nrm   = 0.0;
	norm  = 0.0;

	for (j = 0; j < ncell[1]; ++j )
	{
		y = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			x = xmin[0] + (double(i) + 0.5)*dx[0];
			ij = j * ncell[0] + i;

// Sine function
			err += dx[0]*dx[1]* pow( A[ij] - sin(2.0*pi*x) * sin(2.0*pi*y), 2 );
			nrm += dx[0]*dx[1]* pow( sin(2.0*pi*x) * sin(2.0*pi*y), 2 );

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
		y = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			x = xmin[0] + (double(i) + 0.5)*dx[0];
			ij = j * ncell[0] + i;

			err = A[ij] - sin(2.0*pi*x) * sin(2.0*pi*y);

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

