TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

TRIAD = rgb msg new
DEPEND = $(TRIAD)

include $(TRI)/config/Rules
