TRI = .
ARCH = unix
FFT = fft
POLL = poll

include $(TRI)/config/Common

INCL = 

POLAR = Polar PolarAverage simpfast
RFFT = rfft $(FFT) $(CORE) $(UTILS)
NAVIER = Navier NWave Geometry $(POLAR) Cartesian $(RFFT)
THREEWAVE = ThreeWave NWave Geometry $(POLAR) Cartesian $(RFFT)
BURGER = Burger Geometry Cartesian1 $(RFFT)
POLARAVG = PolarAverageTest PolarAverage simpfast
TRIAD = $(NAVIER) $(THREEWAVE) $(POLARAVG)

burger: $(BURGER:=.o)
	$(C++) $(OPT) -o triad $(BURGER:=.o) $(LIB) 

nw: $(NAVIER:=.o)
	$(C++) $(OPT) -o triad $(NAVIER:=.o) $(LIB) 

w3:	$(THREEWAVE:=.o)
	$(C++) $(OPT) -o triad $(THREEWAVE:=.o) $(LIB)

polaraverage: $(POLARAVG:=.o) $(UTILS:=.o)
	$(C++) $(OPT) -o triad $(POLARAVG=.o) $(UTILS:=.o) $(LIB)

include $(TRI)/config/Rules


