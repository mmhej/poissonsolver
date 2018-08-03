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

	const double c  = 10.0;
	const double r0 = 1.0;
	const double m  = 4.0;

//	const double domain_xmin[3]   = { -2.0, -2.0, -1.0 };
//	const double domain_xmax[3]   = {  2.0,  2.0,  1.0 };

	const double domain_xmin[3]   = { -1.0, -1.0, -1.0 };
	const double domain_xmax[3]   = {  1.0,  1.0,  1.0 };

//	int domain_bounds[3] = { 1, 1, 0 };
	int domain_bounds[3] = { 0, 0, 1 };

	int domain_ncell[3]  = { 64, 64, 32 };
//	int domain_ncell[3]  = { 128, 128, 128 };

//----------------------------------------------------------------------------//
// Variables
//----------------------------------------------------------------------------//
	int     nproc, rank;
	int     d;
	int     i, j, k, kn, kjn, ijk;
	int     n;
	double  x, y, z, r, rho, phi, theta;
	int     ncell[3], icell[3];
	double  dx[3], xmin[3], xmax[3];

	double err, error;
	double nrm, norm;
	double Amag, Bmag;
	double solX, solY, solZ;

	double * Ax;
	double * Ay;
	double * Az;

	double * Bx;
	double * By;
	double * Bz;

	double diffX, diffY, diffZ;

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
// Setup poisson_solver
//----------------------------------------------------------------------------//
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Setup" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	poisson_solver green;
	green.lhs_grad       = true; // specify lhs operator
//	green.lhs_curl       = true; // specify lhs operator
//	green.regularisation = 6;    // regularisation order
//	green.rhs_reproject  = true; // reprojection of rhs field
	green.setup3d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	ncell[0] = green.partition[rank].ncell[0];
	ncell[1] = green.partition[rank].ncell[1];
	ncell[2] = green.partition[rank].ncell[2];

	icell[0] = green.partition[rank].icell[0];
	icell[1] = green.partition[rank].icell[1];
	icell[2] = green.partition[rank].icell[2];

	xmin[0]  = domain_xmin[0] + dx[0] * double(icell[0]);
	xmin[1]  = domain_xmin[1] + dx[1] * double(icell[1]);
	xmin[2]  = domain_xmin[2] + dx[2] * double(icell[2]);

//----------------------------------------------------------------------------//
// Allocate fields
//----------------------------------------------------------------------------//
	Ax    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Ay    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Az    = new double[ncell[0] * ncell[1] * ncell[2]]();

	Bx    = new double[ncell[0] * ncell[1] * ncell[2]]();
	By    = new double[ncell[0] * ncell[1] * ncell[2]]();
	Bz    = new double[ncell[0] * ncell[1] * ncell[2]]();

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

// unbounded-unbounded-periodic
				r = sqrt(x*x + y*y);
				if( r < r0 )
				{
//					Bx[ijk] = pow(1.0 - r*r, m); // Finite order Perlman patch

					Bx[ijk] = 4.0 * c * pow(r0,2) 
					        * exp(- c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
					        * ( pow(r0,4)
					          - pow(x,4)
					          - pow(y,4)
					          - 2.0*pow(x,2)*pow(y,2)
					          - c*pow(x,2)*pow(r0,2)
					          - c*pow(y,2)*pow(r0,2) )
					        * pow(pow(r0,2) - pow(x,2) - pow(y,2),-4);
					By[ijk] = Bx[ijk];
					Bz[ijk] = Bx[ijk];
				}
				else
				{
					Bx[ijk] = 0.0;
					By[ijk] = 0.0;
					Bz[ijk] = 0.0;
				}

// periodic-periodic-unbounded
/*
				if( std::abs(z) < 1.0 )
				{
					Bx[ijk] = ( exp( c*pow(z,2)/((z - 1.0)*(z + 1.0)) )
					          * sin( pi * x)
					          * sin( pi * y)
					          * (  2.0 * pow(pi,2) * pow(z,8)
					            -  8.0 * pow(pi,2) * pow(z,6)
					            + 12.0 * pow(pi,2) * pow(z,4)
					            -  6.0 *         c * pow(z,4)
					            -  8.0 * pow(pi,2) * pow(z,2)
					            -  4.0 *  pow(c,2) * pow(z,2)
					            +  4.0 *         c * pow(z,2)
					            +  2.0 * pow(pi,2)
					            +  2.0 *         c )
					          )*pow( pow(z - 1.0,4) * pow(z + 1.0,4), -1);
					By[ijk] = 0.0;
					Bz[ijk] = 0.0;
				}
				else
				{
					Bx[ijk] = 0.0;
					By[ijk] = 0.0;
					Bz[ijk] = 0.0;
				}
*/
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

	green.push( Bx, By, Bz, NULL, NULL, NULL);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Solve
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Solve" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	green.solve3d( );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Map back to ClientArray
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){ std::cout << " [poisson_solver]: Pull" << std::endl; }
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	green.pull( Bx, By, Bz, Ax, Ay, Az );

//----------------------------------------------------------------------------//
// Calculate error integral
//----------------------------------------------------------------------------//
	err   = 0.0;
	error = 0.0;
	nrm   = 0.0;
	norm  = 0.0;
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

// unbounded-unbounded-periodic
				r = sqrt(x*x + y*y);
				if(r < 0.25*dx[0])
				{
					solX = 0.0;
					solY = 0.0;
					solZ = 0.0;
				}
				else if( r < r0 )
				{
//					solX = - x * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2)); // Perlman patch
//					solY = - y * (1.0 - pow( 1.0 - pow(r,2), m+1) )/(2.0*(m+1.0)*pow(r,2)); // Perlman patch
					solX = - 2.0 * c * pow(r0,2) * x
					     * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
					     * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);
					solY = - 2.0 * c * pow(r0,2) * y
					     * exp( -c * pow(r0,2)/(pow(r0,2) - pow(x,2) - pow(y,2)))
					     * pow(pow(r0,2) - pow(x,2) - pow(y,2), -2);
					solZ = 0.0;
				}
				else
				{
//					solX = - x/(2.0*(m+1.0)*pow(r,2)); // Perlman patch
//					solY = - y/(2.0*(m+1.0)*pow(r,2)); // Perlman patch
					solX = 0.0;
					solY = 0.0;
					solZ = 0.0;
				}

// periodic-periodic-unbounded
/*
				if( std::abs(z) < 1.0 )
				{
					solX = 0.0;
					solY = - 2.0 * c * z * exp( c * pow(z,2)/((z - 1.0)*(z + 1.0)) )
					       * sin( pi * x ) * sin( pi * y )
					       * pow( pow(z - 1.0,2) * pow(z + 1.0,2), -1);
					solZ = - pi * exp( c* pow(z,2)/( (z - 1.0)*(z + 1.0)) )
					       * sin( pi * x ) * cos( pi * y );
				}
				else
				{
					solX = 0.0;
					solY = 0.0;
					solZ = 0.0;
				}
*/

				err += dx[0]*dx[1]*dx[2]* pow( Ax[ijk] - solX, 2);
				nrm += dx[0]*dx[1]*dx[2]* pow( solX, 2);

				err += dx[0]*dx[1]*dx[2]* pow( Ay[ijk] - solY, 2);
				nrm += dx[0]*dx[1]*dx[2]* pow( solY, 2);

				err += dx[0]*dx[1]*dx[2]* pow( Az[ijk] - solZ, 2);
				nrm += dx[0]*dx[1]*dx[2]* pow( solZ, 2);

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
	    << "  " << icell[2] << "  " << icell[2] + ncell[2]
	    <<"' Ghostlevel='0' Origin='"
	    << "  " << domain_xmin[0]
	    << "  " << domain_xmin[1]
	    << "  " << domain_xmin[2]
	    << "' Spacing='"
	    << "  " << dx[0]
	    << "  " << dx[1]
	    << "  " << dx[2] << "'>" << "\n";
	vtifile << "    <Piece Extent='"
	    << "  " << icell[0] << "  " << icell[0] + ncell[0]
	    << "  " << icell[1] << "  " << icell[1] + ncell[1]
	    << "  " << icell[2] << "  " << icell[2] + ncell[2]
	    << "'>" << "\n";
	vtifile << "      <PointData>" << "\n";
	vtifile << "      </PointData>" << "\n";
	vtifile << "      <CellData>" << "\n";

	vtifile << "        <DataArray type='Float64' Name='B' NumberOfComponents='3'  format='ascii'>" << "\n";
	for(ijk = 0; ijk < ncell[0]*ncell[1]*ncell[2]; ++ijk)
	{
		vtifile << std::scientific << std::setw(17) << Bx[ijk] << std::setw(17) << By[ijk] << std::setw(17) << Bz[ijk];
	}
	vtifile << "\n        </DataArray>" << "\n";


	vtifile << "        <DataArray type='Float64' Name='A' NumberOfComponents='3'  format='ascii'>" << "\n";
	for(ijk = 0; ijk < ncell[0]*ncell[1]*ncell[2]; ++ijk)
	{
		vtifile << std::scientific << std::setw(17) << Ax[ijk] << std::setw(17) << Ay[ijk] << std::setw(17) << Az[ijk];
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
		    << "  " << 0 << "  " << domain_ncell[2]
		    <<"' Ghostlevel='0' Origin='"
		    << "  " << domain_xmin[0]
		    << "  " << domain_xmin[1]
		    << "  " << domain_xmin[2]
		    << "' Spacing='"
		    << "  " << dx[0]
		    << "  " << dx[1]
		    << "  " << dx[2] << "'>" << "\n";
		pvtifile << "  <PCellData Vectors='output'>" << "\n";
		pvtifile << "    <PDataArray type='Float64' Name='B' NumberOfComponents='3' format='appended' offset='0'/>" << "\n";
		pvtifile << "    <PDataArray type='Float64' Name='A' NumberOfComponents='3' format='appended' offset='0'/>" << "\n";
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
			         << "  " << icell[2] << "  " << icell[2] + ncell[2]
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

