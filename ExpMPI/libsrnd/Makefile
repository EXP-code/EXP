#
# Extracted libg++ routines
#

ALL: default

CC	= egcc
CCC	= g++
OPT	= -O3 -fPIC    	

OBJ =	ACG.o Binomial.o DiscUnif.o Erlang.o Geom.o HypGeom.o	\
	LogNorm.o MLCG.o NegExp.o Normal.o Poisson.o RNG.o	\
	Random.o Regex.o RndInt.o String.o Uniform.o Weibull.o	\
	error.o timer.o

default: libsrnd.so

libsrnd.so: $(OBJ)
	g++ -shared -o $@ $?
	mv $@ ../


libsrnd.a: $(OBJ)
	ar r $@ $?
	ranlib $@
	mv $@ ../

clean:
	/bin/rm -f *.o *~ *.a *.so

.c.o:
	$(CC) $(OPT) -c $*.c

.cc.o:
	$(CCC) $(OPT) -c $*.cc

