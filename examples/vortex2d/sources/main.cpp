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
#include "class_topology.hpp"
#include "module_output.hpp"

//----------------------------------------------------------------------------//
// MAIN PROGRAM
//----------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
//----------------------------------------------------------------------------//
// Parameters
//----------------------------------------------------------------------------//
	const double pi = acos(-1.0);

	const double r0 = 0.8;
	const double ratio = 2.0; 
	const double c = 2.56085;

	const double domain_xmin[2]   = {-1.0, -1.0};
	const double domain_xmax[2]   = { 1.0,  1.0};

	int domain_bounds[2] = {0,0};

//	int domain_ncell[2]  = { 64, 64 };
	int domain_ncell[2]  = { 128, 128 };

	int     ntime     = 0;
	int     itime     = 0;
	double  dtime     = 0.005;
	double  time      = 0.0;

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

	int stime;

	class_topology topo;

	double * mesh_velx;
	double * mesh_vely;
	double * mesh_vort;

	std::ostringstream str;
	std::string dir, tag;
	std::string filename;
	std::string ifilename;

	std::vector<double> part_posx;
	std::vector<double> part_posy;
	std::vector<double> part_circ;
	std::vector<double> part_velx;
	std::vector<double> part_vely;

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
		std::cout << "EXAMPLE: VORTEX METHOD IN 2D " << std::endl;
		std::cout << std::endl;
		std::cout << "Number of processors intiated: " << nproc << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

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
	poisson.lhs_curl         = true; // specify lhs operator
//	poisson.regularisation = 6;    // regularisation order
	poisson.setup2d( domain_ncell, domain_bounds, dx );

//----------------------------------------------------------------------------//
// Get mesh info
//----------------------------------------------------------------------------//
	topo.ncell[0] = poisson.partition[rank].ncell[0];
	topo.ncell[1] = poisson.partition[rank].ncell[1];

	topo.icell[0] = poisson.partition[rank].icell[0];
	topo.icell[1] = poisson.partition[rank].icell[1];

	topo.xmin[0]  = domain_xmin[0] + dx[0] * double(poisson.partition[rank].icell[0]);
	topo.xmin[1]  = domain_xmin[1] + dx[1] * double(poisson.partition[rank].icell[1]);

	topo.dx[0] = dx[0];
	topo.dx[1] = dx[1];

	topo.domain_ncell[0] = domain_ncell[0];
	topo.domain_ncell[1] = domain_ncell[1];


//----------------------------------------------------------------------------//
// Initiate fields and particles
//----------------------------------------------------------------------------//
	ncell[0] = topo.ncell[0]; ncell[1] = topo.ncell[1];
	dx[0]    = topo.dx[0];    dx[1]    = topo.dx[1];
	xmin[0]  = topo.xmin[0];  xmin[1]  = topo.xmin[1];

	mesh_velx = new double[ncell[0] * ncell[1]]();
	mesh_vely = new double[ncell[0] * ncell[1]]();

	mesh_vort = new double[ncell[0] * ncell[1]]();

	for (j = 0; j < ncell[1]; ++j )
	{
		jn = j * ncell[0];
		y  = xmin[1] + (double(j) + 0.5)*dx[1];
		for (i = 0; i < ncell[0]; ++i )
		{
			ij = jn + i;
			x  = xmin[0] + (double(i) + 0.5)*dx[0];

			r = sqrt( pow(ratio*x,2) + pow(y,2));
			if(r < 0.25*dx[0])
			{
				mesh_vort[ij] = 1.0;

				part_circ.push_back( mesh_vort[ij] * dx[0]*dx[1] );
				part_posx.push_back( x );
				part_posy.push_back( y );

			}
			else if( r < r0 )
			{
				mesh_vort[ij] = ( 1.0 - exp( -c*r0/r * exp( 1.0/(r/r0 - 1.0) ) ) );

				part_circ.push_back( mesh_vort[ij] * dx[0]*dx[1] );
				part_posx.push_back( x );
				part_posy.push_back( y );

			}
			else
			{
				mesh_vort[ij] = 0.0;
			}

		}
	}


	topo.npart = part_circ.size();

//----------------------------------------------------------------------------//
// Simulation loop
//----------------------------------------------------------------------------//
	while( itime <= ntime )
	{

//----------------------------------------------------------------------------//
// Loop over Runge-Kutta steps
//----------------------------------------------------------------------------//
		for( stime = 0; stime < 1; ++stime )
		{


//----------------------------------------------------------------------------//
// Interpolate particles to mesh
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// Solve Poisson equation
//----------------------------------------------------------------------------//
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Push to poisson_solver array
// TO DO: Change to handle ghost cells
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank == 0){ std::cout << " [poisson_solver]: Push" << std::endl; }
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			poisson.push( NULL, NULL, mesh_vort, NULL, NULL, NULL);

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
// TO DO: Change to handle ghost cells
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#ifdef __verb
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank == 0){ std::cout << " [poisson_solver]: Pull" << std::endl; }
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			poisson.pull( NULL, NULL, NULL, mesh_velx, mesh_vely, NULL );


//----------------------------------------------------------------------------//
// Calculate vorticity equation on mesh
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
// Interpolate mesh to particles
//----------------------------------------------------------------------------//




			if(stime == 0)
			{
//----------------------------------------------------------------------------//
// Diagnostics (here to ensure that all values correspond)
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Output field
//----------------------------------------------------------------------------//
#ifdef __verb
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 0){ std::cout << "  Output fields" << std::endl; }
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				str << std::setw(5) << std::setfill('0') << itime;
				dir = "./output";
				tag = "mesh_T" + str.str();
				str.str(""); // clear str

				output_mesh_vtk( dir, tag, topo, mesh_vort );

//----------------------------------------------------------------------------//
// Output particles
//----------------------------------------------------------------------------//
				str << std::setw(5) << std::setfill('0') << itime;
				dir = "./output";
				tag = "part_T" + str.str();
				str.str(""); // clear str

				output_particles( dir, tag, topo, part_posx, part_posy, part_circ );


			} // stime if

//----------------------------------------------------------------------------//
// Advance particles
//----------------------------------------------------------------------------//




		} // RK loop

		itime++;

	} // simulation loop


//----------------------------------------------------------------------------//
// Finalize OpenMPI
//----------------------------------------------------------------------------//
	MPI_Finalize();

	return 0;
}

