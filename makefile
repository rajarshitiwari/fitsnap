
export LAMMPS_DIR=/scratch/softwares/lammps-lib-16Mar18/serial

export MPI_INC=/usr/lib/x86_64-linux-gnu/openmpi/include


LAMMPS-wrapper.o: LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	mpic++ -O3 -I ${MPI_INC} -I ${LAMMPS_DIR} -c LAMMPS-wrapper.cpp

LAMMPS.o: LAMMPS.F90
	mpif90 -O3 -c LAMMPS.F90

single_file.o: single_file.f90
	mpif90 -O3 -c  single_file.f90 -g
fitsnap.x: LAMMPS-wrapper.o LAMMPS.o single_file.o
	mpif90 -I ${MPI_INC} -I ${LAMMPS_DIR} $? -llammps -lmpi_cxx -lstdc++ -lm -g -llapack -o $@ -L ${LAMMPS_DIR}

clean:
	rm *.o *.mod
