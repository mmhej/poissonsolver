//----------------------------------------------------------------------------//
/*
  File:         output_mesh_vtk.cpp

  Description:  Output mesh values to a .vtk file
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// FUNCTION
//----------------------------------------------------------------------------//
void output_mesh_vtk( std::string dir, std::string tag,
                      class_topology topo, double * field )
{
//----------------------------------------------------------------------------//
// Local variables
//----------------------------------------------------------------------------//
	int nproc, rank;
	int i, ij, idx;

	int     domain_ncell[2], ncell[2], icell[2];
	double  dx[2], xmin[2];

	std::ostringstream str;
	std::string filename, ifilename;

//----------------------------------------------------------------------------//
// Get MPI info
//----------------------------------------------------------------------------//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//----------------------------------------------------------------------------//
// Get topology info
//----------------------------------------------------------------------------//
	ncell[0] = topo.ncell[0]; ncell[1] = topo.ncell[1];
	icell[0] = topo.icell[0]; icell[1] = topo.icell[1];
	dx[0]    = topo.dx[0];    dx[1]    = topo.dx[1];
	xmin[0]  = topo.xmin[0];  xmin[1]  = topo.xmin[1];

	domain_ncell[0]  = topo.domain_ncell[0];
	domain_ncell[1]  = topo.domain_ncell[1];

//----------------------------------------------------------------------------//
// Individual .vti files
//----------------------------------------------------------------------------//
	str << std::setw(2) << std::setfill('0') << rank;
	filename = dir + "/" + tag + "_P" + str.str() + ".vti";
	str.str(""); // clear str

	std::ofstream vtifile( filename.c_str() );
//	vtifile.precision(16);
	vtifile.precision(8);
	if(!vtifile.is_open())
	{
		std::cerr << "ERROR: cannot open vtifile." << std::endl;
		return;
	}
	vtifile << "<?xml version='1.0'?>" << "\n";
	vtifile << "<VTKFile type='ImageData' version='0.1' byte_order='LittleEndian'>" << "\n";
	vtifile << "  <ImageData WholeExtent='" 
			<< "  " << icell[0] << "  " << icell[0] + ncell[0]
			<< "  " << icell[1] << "  " << icell[1] + ncell[1]
			<< "  " <<        0 << "  " << 1
			<<"' Ghostlevel='0' Origin='"
			<< "  " << xmin[0]
			<< "  " << xmin[1]
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
	vtifile << "        <DataArray type='Float64' Name='vorticity' NumberOfComponents='1'  format='ascii'>" << "\n";
	for(ij = 0; ij < ncell[0]*ncell[1]; ++ij)
	{
		vtifile << std::scientific << std::setw(17) << field[ij];
	}
	vtifile << "\n        </DataArray>" << "\n";
	vtifile << "      </CellData>" << "\n";
	vtifile << "    </Piece>" << "\n";
	vtifile << "  </ImageData>" << "\n";
	vtifile << "</VTKFile>" << "\n";
	vtifile.close();



//----------------------------------------------------------------------------//
// Main .pvti file
//----------------------------------------------------------------------------//
	filename = dir + "/" + tag + ".pvti";
	std::ofstream pvtifile;

	str << std::setw(2) << std::setfill('0') << rank;
	ifilename = "./" + tag + "_P" + str.str() + ".vti";
	str.str(""); // clear str

	if( rank == 0 )
	{
		pvtifile.open( filename.c_str() );
	//	pvtifile.precision(16);
		pvtifile.precision(8);
		if(!pvtifile.is_open())
		{
			std::cerr << "ERROR: cannot open vtifile." << std::endl;
			return;
		}
		pvtifile << "<?xml version='1.0'?>" << "\n";
		pvtifile << "<VTKFile type='PImageData' version='0.1' byte_order='LittleEndian'>" << "\n";
		pvtifile << "<PImageData WholeExtent='" 
				<< "  " << 0 << "  " << domain_ncell[0]
				<< "  " << 0 << "  " << domain_ncell[1]
				<< "  " << 0 << "  " << 1
				<<"' Ghostlevel='0' Origin='"
				<< "  " << xmin[0]
				<< "  " << xmin[1]
				<< "  " << 0.0
				<< "' Spacing='"
				<< "  " << dx[0]
				<< "  " << dx[1]
				<< "  " << 0.0 << "'>" << "\n";
		pvtifile << "  <PCellData Vectors='output'>" << "\n";
		pvtifile << "    <PDataArray type='Float64' Name='vorticity' NumberOfComponents='1' format='appended' offset='0'/>" << "\n";
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
// Return
//----------------------------------------------------------------------------//
	return;
}

