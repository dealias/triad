TRI = .
ARCH = unix
FFT = fft
POLL = poll

include $(TRI)/config/Common

INCL = 

NAVIER = Navier NWave Polar PolarAverage simpfast Cartesian rfft \
		 $(FFT) $(CORE) $(UTILS)
THREEWAVE = ThreeWave NWave $(CORE) $(UTILS)
BURGER = Burger NWave Cartesian1 rfft $(FFT) $(CORE) $(UTILS)
POLARAVG = PolarAverageTest PolarAverage simpfast

DEPEND = $(NAVIER) $(THREEWAVE) $(BURGER) $(POLARAVG)

nw: $(NAVIER:=.o)
	$(C++) $(OPT) -o triad $(NAVIER:=.o) $(LIB) 

burger: $(BURGER:=.o)
	$(C++) $(OPT) -o triad $(BURGER:=.o) $(LIB) 

w3:	$(THREEWAVE:=.o)
	$(C++) $(OPT) -o triad $(THREEWAVE:=.o) $(LIB)

polaraverage: $(POLARAVG:=.o) $(UTILS:=.o)
	$(C++) $(OPT) -o triad $(POLARAVG=.o) $(UTILS:=.o) $(LIB)

include $(TRI)/config/Rules


