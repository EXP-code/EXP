# 
ALL: default

LAMHOME     = /usr/local/lam62b
SHELL       = /bin/sh
CC          = egcc -D_REENTRANT -I$(LAMHOME)/h
CLINKER     = hf77
PTHREAD     = -lpthread
SYSLIBS     = -L$(LAMHOME)/lib -lmpi $(PTHREAD) -lm
				#Problems using efence with MPI?
# SYSLIBS     = -L$(LAMHOME)/lib -lmpi $(PTHREAD) -lefence -lm
CCC         = g++ -D_REENTRANT  -I$(LAMHOME)/h -I./libsrnd
F77         = g77
CCLINKER    = $(CCC)
OPTFLAGS    = -static -g -DDEBUG
OPTFLAGS    = -O3 -static # -Wall 
CFLAG       = $(OPTFLAGS) -I/home/weinberg/include # -ansi -pedantic
CCFLAG      = $(OPTFLAGS) -I/home/weinberg/include -I/home/weinberg/C++/include -I.
FFLAG       = $(OPTFLAGS) -I/home/weinberg/include

EXECS       = exp


OBJ0 = 		accpot.o begin.o bessacp.o HERNQacp.o CBacp.o end.o	\
		expand.o fixacc.o incpos.o incvel.o init.o initvel.o	\
		outlist.o outdiag.o outlog.o outcoef.o output.o parse.o \
		step.o syncvel.o ddplgndr.o biorth_CB.o biorth_HERNQ.o	\
		biorth_bes.o factrl.o getkey.o potden.o			\
		adjscl.o cubeacp.o ultra.o ultrasphere.o sincos.o	\
		legendre.o tides.o mass_loss.o shock.o pmacp.o		\
	        CBdiskacp.o halo.o outrelx.o pca.o user.o Orient.o	\
		parallel.o Cylacp3.o WghtOrth3.o SatelliteOrbit.o	\
		SLGridMP2.o sledge.o EmpOrth7thd.o biorth1d.o 		\
		SphericalSL.o slabacp.o slabacpSL.o outtipsy.o 		\
		Sphacp.o sincosCPP.o legendreCPP.o parallelCPP.o 	\
		exp_thread.o any_thread2.o scatter.o


default: $(EXECS)

all: default

libsrnd.a:
	cd libsrnd; make OPT="-O3" libsrnd.a

libsrnd.so:
	cd libsrnd; make libsrnd.so

# change to libsrnd.so for shared lib

exp:	$(OBJ0) libsrnd.a
	$(CLINKER) $(OPTFLAGS) -o exp $(OBJ0) \
	-lfft -lnumerical -leigen -lsrnd -lstdc++ -lmdw $(SYSLIBS)

clean:
	/bin/rm -f *.o *~ PI* $(EXECS)
	cd libsrnd; make clean

.c.o:	expand.h
	$(CC) $(CFLAG) -c $*.c

.cc.o:	expand.h
	$(CCC) $(CCFLAG) -c $*.cc

.f.o:
	$(F77) $(FFLAG) -c $*.f


