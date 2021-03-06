export LAMMPS_DIR=/scratch/softwares/lammps-16Mar18

echo ${LAMMPS_DIR}

export MPI_INC=/usr/lib/x86_64-linux-gnu/openmpi/include

mpic++ -O3 -I ${MPI_INC} -I ${LAMMPS_DIR}/src -c LAMMPS-wrapper.cpp 
mpif90 -O3 -c LAMMPS.F90

mpif90 -O3 -c  single_file.f90 -g
mpif90 -I ${MPI_INC} -I ${LAMMPS_DIR}/src LAMMPS.o LAMMPS-wrapper.o single_file.o -llammps_mpi -lmpi_cxx -lstdc++ -lm -g -llapack -o fitsnap.x -L ${LAMMPS_DIR}/src/
