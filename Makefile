LIB = $(LFLAGS) -lm
OPT = -g $(CFLAGS)
MAKEDEPEND = makedepend

ARCH = unix.o fft.o
POLL = poll.o

include config/$(HOSTTYPE)

UTILS = utils.o strcasecmp.o new.o $(POLL) $(ARCH)
CORE = kernel.o Problem.o Integrator.o Param.o $(UTILS)
POLAR = Polar.o PolarAverage.o simpfast.o
NWAVE = NWave.o convolve.o Cartesian.o
TRIAD = Geometry.o $(NWAVE) $(CORE) $(POLAR)


.SUFFIXES: .cc

triad:	Navier.o $(TRIAD)
	$(C++) $(OPT) -o triad Navier.o $(TRIAD) $(LIB) 

w3:	ThreeWave.o $(TRIAD)
	$(C++) $(OPT) -o triad ThreeWave.o $(NWAVE) $(CORE) $(LIB)

kepler:	Kepler.o $(CORE)
	$(C++) $(OPT) -o triad Kepler.o $(CORE) $(LIB)

lotka:	Lotka.o $(CORE)
	$(C++) $(OPT) -o triad Lotka.o $(CORE) $(LIB)

polaraverage: PolarAverageTest.o PolarAverage.o simpfast.o $(UTILS)
	$(C++) $(OPT) -o triad PolarAverageTest.o PolarAverage.o simpfast.o \
		$(UTILS) $(LIB)

.cc.o:
	$(C++) $(OPT) -o $*.o -c $*.cc

clean:
	rm -f *.o *mon.out $(ALL)

depend:
	$(MAKEDEPEND) -f .makedepend $(MDOPT) -I /usr/local/include \
	kernel.cc Problem.cc Integrator.cc Param.cc ThreeWave.cc \
	Navier.cc NWave.cc Geometry.cc Cartesian.cc convolve.cc fft.cc \
	Polar.cc PolarAverage.cc simpfast.cc \
	Kepler.cc Lotka.cc utils.cc strcasecmp.cc new.cc \
	poll.cc idle.cc tremain.cc unix.cc cfft.cc efft.cc

include .makedepend

