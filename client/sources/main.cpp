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

	const double domain_xmin[2]   = {-1.0, -1.0};
	const double domain_xmax[2]   = { 1.0,  1.0};

	int domain_bounds[2] = {0,0};
//	int domain_bounds[2] = {1,1};

//	int domain_ncell[2]  = { 64, 64};
//	int domain_ncell[2]  = { 128, 128};
//	int domain_ncell[2]  = { 1024, 1024};
	int domain_ncell[2]  = { 2048, 2048};
//	int domain_ncell[2]  = { 3000, 3000};
//	int domain_ncell[2]  = { 4096, 4096};

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

	double * test;
	double * sines;
	double * bump;

	double * sines_sol;
	double * bump_sol;

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

//	std::cout << "rank: " << rank << std::endl;
//	MPI_Barrier(MPI_COMM_WORLD);

	for( d = 0; d < 5; ++d)
	{

	domain_ncell[0] = pow(2, d+4);
	domain_ncell[1] = pow(2, d+4);

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
	test  = new double[ncell[0] * ncell[1]]();
	sines = new double[ncell[0] * ncell[1]]();
	bump  = new double[ncell[0] * ncell[1]]();

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

			test[ij] = double(rank);

			sines[ij] = 8.0 * pi * pi * sin(2.0*pi*x) * sin(2.0*pi*y);

			r = sqrt(x*x + y*y);
			if( r < r0 )
			{
				bump[ij] = 4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4);
			}
			else
			{
				bump[ij] = 0.0;
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

	green.push2d( test, sines, bump, NULL, NULL, NULL);

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

	sines_sol = new double[ncell[0] * ncell[1]]();
	bump_sol  = new double[ncell[0] * ncell[1]]();
	green.pull2d( test, sines, bump, NULL, sines_sol, bump_sol);

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
//			err += dx[0]*dx[1]* pow( sines_sol[ij] - sin(2.0*pi*x) * sin(2.0*pi*y), 2 );
//			nrm += dx[0]*dx[1]* pow( sin(2.0*pi*x) * sin(2.0*pi*y), 2 );


// Bump function
			r = sqrt(x*x + y*y);
			if( r < r0 )
			{
				err += dx[0]*dx[1]* pow( bump_sol[ij] - exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2))) , 2);
				nrm += dx[0]*dx[1]* pow( exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2))), 2 );

//vel[1] = - 2.0 * c * pow(r0,2) * y * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);
//vel[2] =   2.0 * c * pow(r0,2) * x * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);

//				err += dx[0]*dx[1]* pow( bump[ij] - (4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4)) ,2);
//				nrm += dx[0]*dx[1]* pow( 4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4) ,2);

			}
			else
			{
				err += dx[0]*dx[1]* pow( bump_sol[ij], 2 );
//				err += dx[0]*dx[1]* pow( bump[ij], 2);
			}


		}
	}

	MPI_Reduce( &err, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce( &nrm, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	if(rank == 0)
	{
		std::cout << "Error: " << std::scientific << std::setw(7) << ncell[0] << std::setw(17) << dx[0] << std::setw(17) << sqrt( error/norm ) << std::endl;
	}


	}



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

//			err = dx[0]*dx[1]* pow( sines_sol[ij] - sin(2.0*pi*x) * sin(2.0*pi*y), 2 );
//			nrm = dx[0]*dx[1]* pow( sin(2.0*pi*x) * sin(2.0*pi*y), 2 );


			r = sqrt(x*x + y*y);
			if( r < r0 )
			{
				err = dx[0]*dx[1]* bump_sol[ij] - exp(-c/(pow(r0,2) - pow(x,2) - pow(y,2)));



//				err = dx[0]*dx[1]* pow( bump[ij] - (4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4)) ,2);
//				nrm = dx[0]*dx[1]* pow( 4.0 * c * pow(r0,2) * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2))) * (pow(r0,4) - pow(x,4) - pow(y,4) - 2.0*pow(x,2)*pow(y,2) - c*pow(x,2)*pow(r0,2) - c*pow(y,2)*pow(r0,2))/pow(pow(r0,2) - pow(x,2) - pow(y,2),4) ,2);

			}
			else
			{
				err = dx[0]*dx[1]* bump_sol[ij];
//				err = dx[0]*dx[1]* pow( bump[ij], 2 );
			}


			outfile << std::scientific << std::setw(17) << x << std::setw(17) << y << std::setw(17) << err << "\n";


//			outfile << std::scientific << std::setw(17) << x << std::setw(17) << y << std::setw(17) << test[ij] << "\n";
//			outfile << std::scientific << std::setw(17) << x << std::setw(17) << y << std::setw(17) << bump[ij] << "\n";

		}
	}
	outfile.close();



//----------------------------------------------------------------------------//
// Finalize OpenMPI
//----------------------------------------------------------------------------//
	MPI_Finalize();

	return 0;
}
