#CC = icc
#CPP = icc
#FC = ifort
#CFLAGS = -m64 -Wall -DGCC_X -O2
#FFLAGS = -m64 -Wall -nofor_main
#LD = icc
#LDFLAGS = -lstdc++

# uncomment for GNU GCC
CC = gcc
CPP = g++
FC = gfortran
LD = g++
CFLAGS = -m64 -Wall -O3
FFLAGS = -m64 -Wall -O3
LDFLAGS = -lgfortran

CPPFLAGS = $(CFLAGS)

#BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLAS = -llapack -lgoto2 -lpthread
BLAS = -framework Accelerate

LDFLAGS += $(BLAS)

MPICC = $(CC) `openmpicc -showme:compile`
MPICPP = $(CPP) `openmpic++ -showme:compile`
MPICFLAGS += -DUSE_MPI
MPIFFLAGS += -DUSE_MPI
MPILD = $(CPP)
MPILDFLAGS = $(LDFLAGS) `openmpic++ -showme:link`

# uncomment on brutus
#MPICC = mpicc
#MPICPP = mpiCC
#MPILD = mpiCC
#MPILDFLAGS = $(LDFLAGS)

