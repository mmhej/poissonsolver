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
  [1] M. M. Hejlesen, G. Winckelmans and J. H. Walther.
      Non-singular Green's functions for the unbounded Poisson equation in one, 
      two and three dimensions.
      arXiv:1704.00704v2 [math.AP], 2018.
  [2] M. M. Hejlesen, J. T. Rasmussen, P. Chatelain, and J. H. Walther. 
      A high order solver for the unbounded Poisson equation. 
      J. Comput. Phys., 252:458--467, 2013.
  [3] H. J. Spietz, M. M. Hejlesen and J. H. Walther.
      A regularization method for solving the Poisson equation for mixed 
      unbounded-periodic domains. 
      J. Comput. Phys., 356:439--447, 2018.


