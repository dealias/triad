TRI = $(HOME)/tri
ARCH = unix
FFT = fft
POLL = poll

include $(TRI)/config/Common

INCL = 

RGB = rgb msg new

DEPEND = $(RGB)

rgb: $(RGB:=.cc)
	$(C++) $(OPT) -o rgb $(RGB:=.cc) $(LIB)

include $(TRI)/config/Rules
