###############################################################################

GPU_OPT=y
DEBUG=n
ADIOS2_OPT=y
MAP_OPT=y
GEM_XGC_COUPLING=y
HOST=crusher

###############################################################################

F90 = ftn

LIBS =/gpfs/alpine/proj-shared/fus123/crusher/GEM_09282022_crusher/dfftpack/libdfftpack.a 

FFTW_DIR =/opt/cray/pe/fftw/3.3.10.1/x86_trento
FFTW_INC =-I$(FFTW_DIR)/include
FFTW_LIB =-L$(FFTW_DIR)/lib -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f

ADIOS2_DIR=/ccs/home/esuchyta/wdmapp/xgc/frontier/update-2022-02-07/adios2-install-cray-8.1.23
ifeq ($(ADIOS2_OPT),y)
	F90 += -D__ADIOS2
	ADIOS2_INC = $(shell $(ADIOS2_DIR)/bin/adios2-config --fortran-flags)
	ADIOS2_LIB = $(shell $(ADIOS2_DIR)/bin/adios2-config --fortran-libs --cxx-libs)
endif

PSPLINE_DIR = /gpfs/alpine/proj-shared/fus123/crusher/GEM_09282022_crusher/pspline
PSPLINE_INC =-I$(PSPLINE_DIR)/include
PSPLINE_LIB =-L$(PSPLINE_DIR)/lib -lpspline

LIBSCI_DIR ?= /opt/cray/pe/libsci/21.08.1.2/amd/40/x86_64
LIBSCI_INC=-I$(LIBSCI_DIR)/include
LIBSCI_LIB=-L$(LIBSCI_DIR)/lib -lsci_amd_mpi -lsci_amd_mpi_mp -lsci_amd -lsci_amd_mp

INC    = $(FFTW_INC) $(ADIOS2_INC)
LD_LIB = $(FFTW_LIB) $(ADIOS2_LIB)


ifeq ($(MAP_OPT),y)
	F90 += -D__MAP_PARALLEL -DEQUIL_OUTPUT
	INC +=$(PSPLINE_INC)
	LD_LIB +=$(PSPLINE_LIB)
ifeq ($(GEM_XGC_COUPLING),y)
	F90 += -D__GEM_XGC_COUPLING
endif
endif

OPT = -O0 -s real64 -hlist=ad -e Zz -homp

ifeq ($(DEBUG), y)
	OPT += -g
endif

ifeq ($(GPU_OPT),y)
	OPT += -hacc -munsafe-fp-atomics -hacc_model=auto_async_none -hacc_model=deep_copy -hacc_model=fast_addr -DGPU
else
	OPT += -DOPENMP_CPU
endif

