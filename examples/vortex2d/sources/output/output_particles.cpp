//----------------------------------------------------------------------------//
/*
  File:         output_particles.cpp

  Description:  Output particle values to file
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// FUNCTION
//----------------------------------------------------------------------------//
void output_particles( std::string dir, std::string tag,
                       class_topology topo, std::vector<double> px,
                       std::vector<double> py, std::vector<double> circ )
{
//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int nproc, rank;
	int i;
	int npart;

	std::ostringstream str;
	std::string filename;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Get topology info
//----------------------------------------------------------------------------//
	npart = topo.npart;

//----------------------------------------------------------------------------//
// Individual .vti files
//----------------------------------------------------------------------------//
	str << std::setw(2) << std::setfill('0') << rank;
	filename = dir + "/" + tag + "_P" + str.str() + ".dat";
	str.str(""); // clear str

	std::ofstream outfile( filename.c_str() );
//	vtifile.precision(16);
	outfile.precision(8);
	if(!outfile.is_open())
	{
		std::cerr << "ERROR: cannot open out file." << std::endl;
		return;
	}
	outfile << "# Particle file" << "\n";
	for(i = 0; i < npart; ++i)
	{
		outfile << std::scientific << std::setw(17) << px[i]
		        << std::scientific << std::setw(17) << py[i]
		        << std::scientific << std::setw(17) << circ[i] << "\n";
	}
	outfile.close();


//----------------------------------------------------------------------------//
// Return
//----------------------------------------------------------------------------//
	return;
}

