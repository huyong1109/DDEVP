#!/bin/tcsh
#PBS -N MPI_TEST
#PBS -q long
#PBS -l nodes=1:ppn=8
#PBS -r n
#PBS -j oe
#PBS -S /bin/csh -V
#cd ${PBS_O_WORKDIR}
source ~/env_machopts.frankfurt_intel
setenv CASE _global
#setenv CASE _DBC
$MPICH_PATH/bin/mpif90 prep$CASE.f90 -L/usr/local/intel-cluster/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95_lp64 -lmkl_blacs_lp64
$MPICH_PATH/bin/mpirun -np 4 ./a.out #>output_prep$CASE
#$MPICH_PATH/bin/mpif77 marching$CASE.f -L/usr/local/intel-cluster/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95_lp64 -lmkl_blacs_lp64
#$MPICH_PATH/bin/mpirun -np 4 ./a.out #>output_marching$CASE

