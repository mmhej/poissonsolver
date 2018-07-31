================================================================================

  GREENFISH
   a Green's function based spectral solver for the Poisson equation for 
   unbounded and/or periodic domains

   Written by Mads Mølholm Hejlesen (2018)

================================================================================


  INTRODUCTION
--------------------------------------------------------------------------------
  GreenFish is a library developed for solving the Poisson equation 
  lap(A) = -B for unbounded and/or periodic domains with spectral accuracy.
  The GreenFish solver is based regularised Green's function solutions which 
  is efficiently computed by an FFT methodology to perform the convolutions.

  The library is distributed under a GNU General Public License by which the 
  terms of usage is specified in the "LICENSE.txt" file included in the main 
  directory.


  INSTALLATION
--------------------------------------------------------------------------------
  The installation of the library is done by a simple make file in the main 
  directory. If the make command does not work directly you may have to manually
  edit the "makefile" file.
  The only installation requirements at this point is OpenMPI.

  NB. The library is developed on a LINUX platform and problems may arise if 
      used on another platform.


  LANGUAGE INTERFACES
--------------------------------------------------------------------------------
  At the present moment the only interface to the GreenFish library is in C++.
  Future plans include interfaces to C, Fortran and Python.


  EXAMPLES AND TUTORIALS
--------------------------------------------------------------------------------
  In the "examples" directory a number of examples are included which may serve 
  for now serve as the main tutorial.


  LITERATURE
--------------------------------------------------------------------------------
  



