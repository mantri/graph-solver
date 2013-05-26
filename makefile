MKLROOT=/opt/apps/intel/mkl/10.2.2.025
INTELROOT=/opt/apps/intel/11.1/

#CLIBS=-Wl,--start-group -L$(MKLROOT)/lib/em64t -L$(INTELROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -Wl,--end-group -limf -lsvml -lirc -liomp5 -lpthread -lm
CLIBS=-Wl,--start-group -L$(MKLROOT)/lib/em64t -L$(INTELROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -Wl,--end-group  -lpthread -lm -lrt -liomp5

LIB_LAPACK = -L$(MKLROOT)/lib/64 -L$(TACC_MKL_LIB) -L/opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/

LIB = $(CLIBS) $(LIB_LAPACK)

#CFLAGS = -m64 -O0  -g -openmp -std=c99 -pthread
CFLAGS = -m64 -O0  -g -openmp 
INCL_MKL =  -I$(MKLROOT)/include -I$(TACC_MKL_INC)

all: 
	icc $(CFLAGS) -c modifiedGS.c
	icc -o MGS modifiedGS.o $(LIB)
	icc $(CFLAGS) $(INCL_MKL) -c mlgps.c
	icc -o MLGPS mlgps.o $(LIB) 
	icc $(CFLAGS) $(INCL_MKL) -c mlgps_mlcb.c
	icc -o MLGPS_MLCB mlgps_mlcb.o $(LIB) 
	rm -rf *.o

clean:
	rm MLGPS MLGPS_MLCB MGS
