

# Installing fftw with enable-mpi option

# Download and unpack

./configure CC=gcc CXX=g++ --prefix=/home/mmhej/Programs/fftw-3.3.5_gcc-4.8 --enable-mpi

make

make check

make install

make clean
