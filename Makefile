TRI = .
ARCH = unix
FFT = fft rfft
POLL = poll

include $(TRI)/config/Common

INCL = 

POLAR = Polar PolarAverage simpfast
NWAVE = NWave Geometry convolve Cartesian $(FFT) $(CORE) $(UTILS)
NAVIER = Navier $(NWAVE) $(POLAR)
THREEWAVE = ThreeWave $(NWAVE)
POLARAVG = PolarAverageTest PolarAverage simpfast
TRIAD = $(NAVIER) $(THREEWAVE) $(POLARAVG)

triad: $(NAVIER:=.o)
	$(C++) $(OPT) -o triad $(NAVIER:=.o) $(LIB) 

w3:	$(THREEWAVE:=.o)
	$(C++) $(OPT) -o triad $(THREEWAVE:=.o) $(LIB)

polaraverage: $(POLARAVG:=.o) $(UTILS:=.o)
	$(C++) $(OPT) -o triad $(POLARAVG=.o) $(UTILS:=.o) $(LIB)

include $(TRI)/config/Rules


