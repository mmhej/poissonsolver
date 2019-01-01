//-------------------------------------------------------------------------------//
/*
  File:         module_mesh.hpp

  Description:  Header file for module_mesh.cpp
*/
//-------------------------------------------------------------------------------//
#ifndef __incl_output
#define __incl_output

//-------------------------------------------------------------------------------//
// HEADERS INCLUDE
//-------------------------------------------------------------------------------//
// System
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
// External
#include "mpi.h"
// Internal
#include "class_topology.hpp"

//-------------------------------------------------------------------------------//
// PROTOTYPE FUNCTIONS
//-------------------------------------------------------------------------------//
void output_mesh_vtk( std::string, std::string, class_topology, double * );

void output_particles( std::string, std::string, class_topology,
                       std::vector<double>, std::vector<double>, std::vector<double> );

#endif
