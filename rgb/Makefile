TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

RGB = rgb palette msg new getopt
DEPEND = $(RGB)

rgb: $(RGB:=.o)
	$(C++) $(OPT) -o rgb $(RGB:=.o) $(LIB)

include $(TRI)/config/Rules
