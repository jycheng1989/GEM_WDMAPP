#OPTION
GPU_OPT=y
DEBUG=n
ADIOS2_OPT=y
MAP_OPT=y
GEM_XGC_COUPLING=y
HOST=crusher
###############################################################################
ifneq (,${HOST})
  SYSTEMS := ${HOST}
else
  SYSTEMS := $(shell hostname)
endif

ifneq (,$(findstring cori,$(SYSTEMS)))
	GPU_OPT=n
endif

SRCS := $(wildcard *.F90)
OBJS := $(patsubst %.F90,%.o,$(SRCS))
OBJS =  gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o
ifeq ($(ADIOS2_OPT),y)
	SRCS += adios2_comm_mod.F90
endif
ifeq ($(ADIOS2_OPT),y)
	OBJS += adios2_comm_mod.o
endif
ifeq ($(MAP_OPT),y)
	OBJS += mapping.o
endif
ifeq ($(GEM_XGC_COUPLING),y)
	OBJS += coupling_core_edge.o
endif
ifneq (,$(findstring summit,$(SYSTEMS)))
   LIBS = /ccs/home/cycheng/software/nvhpc21.7/dfftpack/libdfftpack.a 
   FFTW_DIR ?=/sw/summit/spack-envs/base/opt/linux-rhel8-ppc64le/nvhpc-21.7/fftw-3.3.9-bzi7deue27ijd7xm4zn7pt22u4sj47g4
endif

ifneq (,$(findstring crusher,$(SYSTEMS)))
  LIBS =/ccs/home/cycheng/software/crusher/cray/dfftpack/libdfftpack.a 
  FFTW_DIR ?=/opt/cray/pe/fftw/3.3.8.13/x86_trento
endif

ifneq (,$(findstring cori,$(SYSTEMS)))
	LIBS = /project/projectdirs/mp118/jycheng/SOFTWARE/dfftpack-mic/libdfftpack.a
	FFTW_DIR = /opt/cray/pe/fftw/default/mic_knl
endif

FFTW_INC ?= -I$(FFTW_DIR)/include
FFTW_LIB ?= -L$(FFTW_DIR)/lib -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f

#ifneq (,$(findstring cori,$(SYSTEMS)))
#	FFTW_INC ?= -I/opt/cray/pe/fftw/default/mic_knl/include -I/opt/cray/pe/fftw/default/mic_knl/mod
#	FFTW_LIB ?= -L/opt/cray/pe/fftw/default/mic_knl/lib -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f
#endif

ifeq ($(ADIOS2_OPT),y)
ifneq (,$(findstring summit,$(SYSTEMS)))
	#ADIOS2_DIR ?= /gpfs/alpine/world-shared/phy122/lib/install/summit/adios2/devel/nvhpc
    #ADIOS2_DIR ?= /gpfs/alpine/world-shared/fus123/jycheng/SOFTWARE/nvhpc21.7/ADIOS2-bin
    ADIOS2_DIR ?=  /ccs/home/merlo/soft/adios2_XGC_RDMA
	#ADIOS2_DIR ?= /gpfs/alpine/world-shared/phy122/lib/install/summit/adios2/devel/nvhpc
endif

ifneq (,$(findstring crusher,$(SYSTEMS)))
    #ADIOS2_DIR ?= /autofs/nccs-svm1_sw/crusher/spack-envs/base/opt/cray-sles15-zen3/cce-14.0.2/adios2-2.8.1-5qitwvfyf4m3jtzxq6uaczvjk6prhwkj
    #ADIOS2_DIR ?= /autofs/nccs-svm1_sw/crusher/spack-envs/base/opt/cray-sles15-zen3/cce-14.0.3/adios2-2.8.1-w6s7ykstn5qpq66rggkueewlmdkqclap
    ADIOS2_DIR ?= /gpfs/alpine/proj-shared/fus123/jycheng/merge-07052022/GEM_09282022_crusher/ADIOS2/adios2_install
endif

ifneq (,$(findstring cori,$(SYSTEMS)))
	ADIOS2_DIR ?= /project/projectdirs/m499/Software/adios2/DEFAULT/cori_knl/intel-static
endif
	ADIOS2_INC ?= $(shell $(ADIOS2_DIR)/bin/adios2-config --fortran-flags)
	ADIOS2_LIB ?= $(shell $(ADIOS2_DIR)/bin/adios2-config --fortran-libs --cxx-libs)
endif

ifeq ($(MAP_OPT),y)
ifneq (,$(findstring summit,$(SYSTEMS)))
	PSPLINE_DIR ?= /gpfs/alpine/world-shared/phy122/lib/install/summit/pspline/nvhpc21.7
endif
ifneq (,$(findstring cori,$(SYSTEMS)))
	PSPLINE_DIR ?= /project/projectdirs/m499/Software/pspline/DEFAULT/cori_knl/DEFAULT
endif
ifneq (,$(findstring crusher,$(SYSTEMS)))
    PSPLINE_DIR ?= /ccs/home/cycheng/software/crusher/cray/pspline
endif
	PSPLINE_INC=-I$(PSPLINE_DIR)/include
	PSPLINE_LIB=-L$(PSPLINE_DIR)/lib -lpspline
endif

ifneq (,$(findstring crusher,$(SYSTEMS)))
    LIBSCI_DIR ?= /opt/cray/pe/libsci/21.08.1.2/amd/40/x86_64
    LIBSCI_INC=-I$(LIBSCI_DIR)/include
    LIBSCI_LIB=-L$(LIBSCI_DIR)/lib -lsci_amd_mpi -lsci_amd_mpi_mp -lsci_amd -lsci_amd_mp
endif

LIB =
LD_LIB =

ifeq ($(ADIOS2_OPT),y)
	LIB +=$(ADIOS2_INC) #-I/autofs/nccs-svm1_sw/summit/cuda/10.2.89/include -I/autofs/nccs-svm1_sw/summit/cuda/10.2.89/mod
	LD_LIB +=$(ADIOS2_LIB) #-L/autofs/nccs-svm1_sw/summit/cuda/10.2.89/lib64 -lcublas
endif

ifneq (,$(findstring cori,$(SYSTEMS)))
LIB +=-I/opt/cray/pe/fftw/default/mic_knl/include -I/opt/cray/pe/fftw/default/mic_knl/mod
LD_LIB +=-L/opt/cray/pe/fftw/default/mic_knl/lib -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f
endif

ifneq (,$(findstring summit,$(SYSTEMS)))
LIB += $(FFTW_INC)
LD_LIB += $(FFTW_LIB)
endif

ifneq (,$(findstring crusher,$(SYSTEMS)))
LIB +=$(FFTW_INC)
LD_LIB +=$(FFTW_LIB)
#LIB +=-I/opt/rocm-5.1.0/rocblas/include
#LIB +=-I/opt/cray/pe/libsci/default/AMD/40/x86_64/include
#LD_LIB +=-L/opt/rocm-5.1.0/rocblas/lib -lrocblas
#LD_LIB +=-L/opt/cray/pe/libsci/default/AMD/40/x86_64/lib/libsci_amd.a
#LIB +=$(LIBSCI_INC)
#LD_LIB +=$(LIBSCI_LIB)
endif


ifneq (,$(findstring summit,$(SYSTEMS)))
F90 = mpif90
#export TAU_MAKEFILE=/ccs/home/cycheng/software/nvhpc21.7/Tau/sourcecode/tau-openmp-intall/tau-install/ibm64linux/lib/Makefile.tau-acc-mpi-openmp-nvhpc
#F90 = /ccs/home/cycheng/software/nvhpc21.7/Tau/sourcecode/tau-openmp-intall/tau-install/ibm64linux/bin/tau_f90.sh -tau_options=-optCompInst
endif
ifneq (,$(findstring cori,$(SYSTEMS)))
F90 = ftn
endif
ifneq (,$(findstring crusher,$(SYSTEMS)))
F90 = ftn
endif
PLIB = gem_pputil.o

ifeq ($(ADIOS2_OPT),y)
	F90 += -D__ADIOS2
endif

ifeq ($(MAP_OPT),y)
	F90 += -D__MAP_PARALLEL -DEQUIL_OUTPUT
	LIB +=$(PSPLINE_INC)
	LD_LIB +=$(PSPLINE_LIB)
ifeq ($(GEM_XGC_COUPLING),y)
	F90 += -D__GEM_XGC_COUPLING
endif
endif
#OPT = -f free -s real64 -eD -h omp -Ktrap=fp -m 4  
#OPT = -f free -s real64 -eD -h omp
#OPT = -f free -s real64 -O2 -h omp

ifneq (,$(findstring summit,$(SYSTEMS)))
OPT = -O0 -r8 -Kieee -llapack -lblas -g -cpp -Minfo=accel -acc -ta=nvidia:cc70
ifeq ($(GPU_OPT),y)
	OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Minfo=accel -acc -mp -ta=nvidia:cc70 #-Mcuda -Mcudalib=cublas #-hlist=a #-lcudart -lcublas 
ifeq ($(DEBUG), y)
	OPT += -g -Mbounds
endif
	OPT += -DGPU 
else
	OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -mp
ifeq ($(DEBUG),y)
	OPT += -g -Mbounds
endif
	OPT += -DOPENMP_CPU
endif
endif

ifneq (,$(findstring crusher,$(SYSTEMS)))
OPT = -O0-h zero -s real64 -hlist=ad -e Zz -homp
ifeq ($(GPU_OPT),y)
    OPT = -O0 -h zero -s real64 -hlist=ad -e Zz -hacc -munsafe-fp-atomics -hacc_model=auto_async_none
ifeq ($(DEBUG), y)
    OPT += -g 
endif
    OPT += -DGPU
else
    OPT = -O0 -h zero -s real64 -hlist=ad -e Zz -homp 
ifeq ($(DEBUG),y)
    OPT += -g
endif
    OPT += -DOPENMP_CPU
endif
endif

ifneq (,$(findstring cori,$(SYSTEMS)))
	OPT = -O3 -FR -r8 -qopenmp
ifeq ($(DEBUG), y)
	OPT += -g -traceback -check bounds
endif
	OPT += -DOPENMP_CPU
endif

#all : gem
ifeq ($(ADIOS2_OPT),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(LIB) $(LD_LIB) 
else
ifeq ($(MAP_OPT),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(LIB) $(LD_LIB)
else
ifeq ($(GEM_XGC_COUPLING),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o mapping.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(LIB) $(LD_LIB)
else
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o mapping.o coupling_core_edge.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(LIB) $(LD_LIB)
endif
endif
endif

gem_pputil.o: gem_pputil.F90
	$(F90) -c $(OPT) gem_pputil.F90

gem_com.o: gem_com.F90 gem_pputil.o
	$(F90) -c $(OPT) gem_com.F90

gem_equil.o: gem_equil.F90 gem_pputil.o gem_com.o
	$(F90) -c $(OPT) gem_equil.F90

gem_gkps_adi.o: gem_gkps_adi.F90 gem_com.F90 gem_equil.F90 gem_pputil.F90
	$(F90) -c $(OPT) gem_gkps_adi.F90

ifeq ($(ADIOS2_OPT),y)
adios2_comm_mod.o : adios2_comm_mod.F90 gem_com.o
	$(F90) -c $(OPT) $(LIB) adios2_comm_mod.F90
endif

ifeq ($(MAP_OPT),y)
mapping.o : mapping.F90 gem_com.o adios2_comm_mod.o
	$(F90) -c $(OPT) $(LIB) mapping.F90
ifeq ($(GEM_XGC_COUPLING),y)
coupling_core_edge.o :coupling_core_edge.F90 gem_com.o gem_equil.o adios2_comm_mod.o mapping.o
	$(F90) -c $(OPT) $(LIB) coupling_core_edge.F90
endif
endif

ifeq ($(ADIOS2_OPT),n)
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o
	$(F90) -c $(OPT) $(LIB) gem_main.F90
else
ifeq ($(MAP_OPT),n)
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o
	$(F90) -c $(OPT) $(LIB) gem_main.F90
else
ifeq ($(GEM_XGC_COUPLING),n)
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o mapping.o
	$(F90) -c $(OPT) $(LIB) gem_main.F90
else
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o mapping.o coupling_core_edge.o
	$(F90) -c $(OPT) $(LIB) gem_main.F90
endif
endif
endif

gem_outd.o: gem_outd.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) gem_outd.F90

gem_fcnt.o: gem_fcnt.F90
	$(F90) -c $(OPT) gem_fcnt.F90

gem_fft_wrapper.o: gem_fft_wrapper.F90
	$(F90) -c $(OPT) $(LIB) gem_fft_wrapper.F90

clean:
	rm -f *.o *.opt *.i *acc.o *acc.s *.lst *.mod *.cg gem_main
