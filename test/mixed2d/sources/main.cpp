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
#include "poisson_solver.hpp"
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

	const double domain_xmin[2]   = {-1.0, -2.0};
	const double domain_xmax[2]   = { 1.0,  2.0};

	int domain_bounds[2] = {1,0};

//	int domain_ncell[2]  = { 64, 64};
//	int domain_ncell[2]  = { 128, 128};

	int domain_ncell[2]  = { 128, 256};

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
	int     ncell[2], icell[2];
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
	std::string ifilename;

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
		std::cout << "TESTING MIXED PERIODIC-UNBOUNDED DOMAIN IN 2D... " << std::endl;
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
// Setup poisson_solver
//----------------------------------------------------------------------------//
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Setup" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	poisson_solver poisson;
	poisson.lhs_grad       = true; // specify lhs operator
//	poisson.regularisation = 6;    // regularisation order
	poisson.setup2d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	ncell[0] = poisson.partition[rank].ncell[0];
	ncell[1] = poisson.partition[rank].ncell[1];

	icell[0] = poisson.partition[rank].icell[0];
	icell[1] = poisson.partition[rank].icell[1];

	xmin[0]  = domain_xmin[0] + dx[0] * double(poisson.partition[rank].icell[0]);
	xmin[1]  = domain_xmin[1] + dx[1] * double(poisson.partition[rank].icell[1]);

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
		y = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x = xmin[0] + (double(i) + 0.5)*dx[0];

			if(std::abs(y) < r0)
			{
/*
				Bx[ij] = pow(1.0 - y*y, 5);
*/


				Bx[ij] = ( exp( c * pow(y,2) / ( (y - 1.0)*(y + 1.0) ) )
				       * sin( pi * x ) 
				       * ( 1.0 * pow(pi,2) * pow(y,8) 
				         - 4.0 * pow(pi,2) * pow(y,6) 
				         + 6.0 * pow(pi,2) * pow(y,4)
				         - 6.0 * c         * pow(y,4) 
				         - 4.0 * pow(pi,2) * pow(y,2)
				         - 4.0 * pow(c,2)  * pow(y,2)
				         + 4.0 * c         * pow(y,2) 
				         + 1.0 * pow(pi,2) 
				         + 2.0 * c ) 
				         )/( pow(y - 1.0,4) * pow(y + 1.0,4) );

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
// Push to poisson_solver array
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Push" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	poisson.push( Bx, NULL, NULL, NULL, NULL, NULL);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Solve
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Solve" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	poisson.solve2d( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Map back to ClientArray
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Pull" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	poisson.pull( Bx, NULL, NULL, Ax, Ay, NULL );

//----------------------------------------------------------------------------//
// Calculate error
//----------------------------------------------------------------------------//
	err   = 0.0;
	error = 0.0;
	nrm   = 0.0;
	norm  = 0.0;

	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
		y = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x = xmin[0] + (double(i) + 0.5)*dx[0];

			if(std::abs(y) < r0)
			{
/*
				solX = 0.0;
				solY = 1.0/11.0*pow(y,11) - 5.0/9.0*pow(y,9) + 10.0/7.0*pow(y,7) - 2.0*pow(y,5) + 5.0/3.0*pow(y,3) - y;
*/

				solX = pi * exp( c * pow(y,2)/( ( y - 1.0 )*( y + 1.0 ) ) ) * cos( pi*x );
				solY = - 2.0 * c * y * exp( c * pow(y,2) /( ( y - 1.0 )*( y + 1.0 ) ) ) * sin( pi * x )/( pow( y - 1.0, 2 ) * pow( y + 1.0, 2 ) );

			}
			else
			{
/*
				if( y > 0.0)
				{
					solX = 0.0;
					solY = - 256.0/693.0;
				}
				else if( y < 0.0)
				{
					solX = 0.0;
					solY = 256.0/693.0;
				}
*/

				solX = 0.0;
				solY = 0.0;

			}

			err += dx[0]*dx[1]* pow( Ax[ij] - solX , 2);
			nrm += dx[0]*dx[1]* pow( solX, 2);

			err += dx[0]*dx[1]* pow( Ay[ij] - solY , 2);
			nrm += dx[0]*dx[1]* pow( solY, 2);


		}
	}

	MPI_Reduce( &err, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce( &nrm,  &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Individual .vti files
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	str << std::setw(2) << std::setfill('0') << rank;
	filename = "./output/mesh_P" + str.str() + ".vti";
	str.str(""); // clear str
	std::ofstream vtifile( filename.c_str() );
//	vtifile.precision(16);
	vtifile.precision(8);
	if(!vtifile.is_open())
	{
		std::cerr << "ERROR: cannot open vtifile." << std::endl;
		return 0;
	}
	vtifile << "<?xml version='1.0'?>" << "\n";
	vtifile << "<VTKFile type='ImageData' version='0.1' byte_order='LittleEndian'>" << "\n";
	vtifile << "  <ImageData WholeExtent='" 
	    << "  " << icell[0] << "  " << icell[0] + ncell[0]
	    << "  " << icell[1] << "  " << icell[1] + ncell[1]
	    << "  " <<        0 << "  " << 1
	    <<"' Ghostlevel='0' Origin='"
	    << "  " << domain_xmin[0]
	    << "  " << domain_xmin[1]
	    << "  " << 0.0
	    << "' Spacing='"
	    << "  " << dx[0]
	    << "  " << dx[1]
	    << "  " << 0.0 << "'>" << "\n";
	vtifile << "    <Piece Extent='"
	    << "  " << icell[0] << "  " << icell[0] + ncell[0]
	    << "  " << icell[1] << "  " << icell[1] + ncell[1]
	    << "  " <<        0 << "  " << 1
	    << "'>" << "\n";
	vtifile << "      <PointData>" << "\n";
	vtifile << "      </PointData>" << "\n";
	vtifile << "      <CellData>" << "\n";
	vtifile << "        <DataArray type='Float64' Name='B' NumberOfComponents='1'  format='ascii'>" << "\n";
	for(ij = 0; ij < ncell[0]*ncell[1]; ++ij)
	{
		vtifile << std::scientific << std::setw(17) << Bx[ij];
	}
	vtifile << "\n        </DataArray>" << "\n";
	vtifile << "        <DataArray type='Float64' Name='A' NumberOfComponents='2'  format='ascii'>" << "\n";
	for(ij = 0; ij < ncell[0]*ncell[1]; ++ij)
	{
		vtifile << std::scientific << std::setw(17) << Ax[ij] << std::setw(17) << Ay[ij];
	}
	vtifile << "\n        </DataArray>" << "\n";
	vtifile << "      </CellData>" << "\n";
	vtifile << "    </Piece>" << "\n";
	vtifile << "  </ImageData>" << "\n";
	vtifile << "</VTKFile>" << "\n";
	vtifile.close();

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Main .pvti file
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	filename = "./output/mesh.pvti";
	std::ofstream pvtifile;

	str << std::setw(2) << std::setfill('0') << rank;
	ifilename = "./mesh_P" + str.str() + ".vti";
	str.str(""); // clear str

	if( rank == 0 )
	{
		pvtifile.open( filename.c_str() );
	//	pvtifile.precision(16);
		pvtifile.precision(8);
		if(!pvtifile.is_open())
		{
			std::cerr << "ERROR: cannot open vtifile." << std::endl;
			return 0;
		}
		pvtifile << "<?xml version='1.0'?>" << "\n";
		pvtifile << "<VTKFile type='PImageData' version='0.1' byte_order='LittleEndian'>" << "\n";
		pvtifile << "<PImageData WholeExtent='" 
		    << "  " << 0 << "  " << domain_ncell[0]
		    << "  " << 0 << "  " << domain_ncell[1]
		    << "  " << 0 << "  " << 1
		    <<"' Ghostlevel='0' Origin='"
		    << "  " << domain_xmin[0]
		    << "  " << domain_xmin[1]
		    << "  " << 0.0
		    << "' Spacing='"
		    << "  " << dx[0]
		    << "  " << dx[1]
		    << "  " << 0.0 << "'>" << "\n";
		pvtifile << "  <PCellData Vectors='output'>" << "\n";
		pvtifile << "    <PDataArray type='Float64' Name='B' NumberOfComponents='1' format='appended' offset='0'/>" << "\n";
		pvtifile << "    <PDataArray type='Float64' Name='A' NumberOfComponents='2' format='appended' offset='0'/>" << "\n";
		pvtifile << "  </PCellData>" << "\n";
		pvtifile.close();
	}

	for(i = 0; i < nproc; ++i)
	{
		if( rank == i )
		{
			pvtifile.open( filename.c_str(), std::ofstream::app );
			pvtifile << "  <Piece Extent='"
			         << "  " << icell[0] << "  " << icell[0] + ncell[0]
			         << "  " << icell[1] << "  " << icell[1] + ncell[1]
			         << "  " <<        0 << "  " << 1
			         << "' Source='" << ifilename << "'/>" << "\n";
			if(rank == nproc - 1)
			{
				pvtifile << "</PImageData>" << "\n";
				pvtifile << "</VTKFile>" << "\n";
			}
			pvtifile.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

//----------------------------------------------------------------------------//
// Finalize OpenMPI
//----------------------------------------------------------------------------//
	MPI_Finalize();

	return 0;
}

