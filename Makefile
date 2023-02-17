HOST=crusher
include $(HOST).mk

DEBUG ?= n
ADIOS2_OPT ?= y
GPU_OPT ?= y
MAP_OPT ?= y
GEM_XGC_COUPLING ?= y

F90 ?= mpif90


SRCS := $(wildcard *.F90)
OBJS := $(patsubst %.F90,%.o,$(SRCS))
OBJS =  gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o
ifeq ($(MAP_OPT),y)
	OBJS += mapping.o
endif
ifeq ($(GEM_XGC_COUPLING),y)
	OBJS += coupling_core_edge.o
endif
PLIB = gem_pputil.o


#all : gem
ifeq ($(ADIOS2_OPT),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(INC) $(LD_LIB) 
else
ifeq ($(MAP_OPT),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(INC) $(LD_LIB)
else
ifeq ($(GEM_XGC_COUPLING),n)
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o mapping.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(INC) $(LD_LIB)
else
gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o adios2_comm_mod.o mapping.o coupling_core_edge.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) $(INC) $(LD_LIB)
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

adios2_comm_mod.o : adios2_comm_mod.F90 gem_com.o
	$(F90) -c $(OPT) $(INC) adios2_comm_mod.F90

ifeq ($(MAP_OPT),y)
mapping.o : mapping.F90 gem_com.o adios2_comm_mod.o
	$(F90) -c $(OPT) $(INC) mapping.F90
ifeq ($(GEM_XGC_COUPLING),y)
coupling_core_edge.o :coupling_core_edge.F90 gem_com.o gem_equil.o adios2_comm_mod.o mapping.o
	$(F90) -c $(OPT) $(INC) coupling_core_edge.F90
endif
endif

ifeq ($(MAP_OPT),n)
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o mapping.o
	$(F90) -c $(OPT) $(INC) gem_main.F90
else
ifeq ($(GEM_XGC_COUPLING),n)
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o mapping.o
	$(F90) -c $(OPT) $(INC) gem_main.F90
else
gem_main.o: gem_main.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o adios2_comm_mod.o mapping.o coupling_core_edge.o
	$(F90) -c $(OPT) $(INC) gem_main.F90
endif
endif

gem_outd.o: gem_outd.F90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) gem_outd.F90

gem_fcnt.o: gem_fcnt.F90
	$(F90) -c $(OPT) gem_fcnt.F90

gem_fft_wrapper.o: gem_fft_wrapper.F90
	$(F90) -c $(OPT) $(INC) gem_fft_wrapper.F90

clean:
	rm -f *.o *.opt *.i *acc.o *acc.s *.lst *.mod *.cg gem_main
