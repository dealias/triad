G++ = g++ -Wall
LIB = $(LFLAGS) -lm
OPT = -g $(CFLAGS)
MAKEDEPEND = makedepend

POLL = poll.o

include config/$(HOSTTYPE)

UTILS = utils.o strcasecmp.o new.o $(POLL)
CORE = kernel.o Approx.o Integrator.o Param.o $(UTILS)
POLAR = Polar.o PolarAverage.o simpfast.o
TRIAD = $(CORE) NWave.o Geometry.o $(POLAR)

.SUFFIXES: .cc

triad:	Navier.o $(TRIAD)
	$(G++) $(OPT) -o triad Navier.o $(TRIAD) $(LIB) 

w3:	ThreeWave.o $(TRIAD)
	$(G++) $(OPT) -o triad ThreeWave.o NWave.o $(CORE) $(LIB)

kepler:	Kepler.o $(CORE)
	$(G++) $(OPT) -o triad Kepler.o $(CORE) $(LIB)

lotka:	Lotka.o $(CORE)
	$(G++) $(OPT) -o triad Lotka.o $(CORE) $(LIB)

polaraverage: PolarAverageTest.o PolarAverage.o simpfast.o $(UTILS)
	$(G++) $(OPT) -o triad PolarAverageTest.o PolarAverage.o simpfast.o \
		$(UTILS) $(LIB)

.cc.o:
	$(G++) $(OPT) -o $*.o -c $*.cc

clean:
	rm -f *.o *mon.out $(ALL)

depend:
	$(MAKEDEPEND) $(MDOPT) -I /usr/local/include \
	kernel.cc Approx.cc Integrator.cc Param.cc ThreeWave.cc \
	Navier.cc NWave.cc Geometry.cc Polar.cc PolarAverage.cc simpfast.cc \
	Kepler.cc Lotka.cc utils.cc strcasecmp.cc new.cc idle.cc

# DO NOT DELETE THIS LINE -- make depend depends on it.

kernel.o: /usr/lpp/xlC/include/iomanip.h /usr/lpp/xlC/include/generic.h
kernel.o: /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
kernel.o: /usr/include/string.h /usr/include/standards.h
kernel.o: /usr/include/sys/types.h kernel.h /usr/lpp/xlC/include/fstream.h
kernel.o: /usr/include/stdio.h /usr/include/limits.h
kernel.o: /usr/include/sys/limits.h /usr/lpp/xlC/include/math.h
kernel.o: /usr/include/math.h /usr/include/errno.h /usr/include/time.h
kernel.o: types.h precision.h /usr/include/float.h Complex.h utils.h
kernel.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
kernel.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
kernel.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
kernel.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
kernel.o: Table.h Param.h Integrator.h Approx.h
Approx.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Approx.o: /usr/include/string.h /usr/include/standards.h
Approx.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Approx.o: /usr/include/stdio.h /usr/include/limits.h
Approx.o: /usr/include/sys/limits.h /usr/lpp/xlC/include/math.h
Approx.o: /usr/include/math.h /usr/include/errno.h /usr/include/time.h
Approx.o: types.h precision.h /usr/include/float.h Complex.h utils.h
Approx.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
Approx.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
Approx.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
Approx.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
Approx.o: Table.h Param.h Integrator.h Approx.h
Integrator.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Integrator.o: /usr/include/string.h /usr/include/standards.h
Integrator.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Integrator.o: /usr/include/stdio.h /usr/include/limits.h
Integrator.o: /usr/include/sys/limits.h /usr/lpp/xlC/include/math.h
Integrator.o: /usr/include/math.h /usr/include/errno.h /usr/include/time.h
Integrator.o: types.h precision.h /usr/include/float.h Complex.h utils.h
Integrator.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
Integrator.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
Integrator.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
Integrator.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
Integrator.o: Table.h Param.h Integrator.h Approx.h
Param.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Param.o: /usr/include/string.h /usr/include/standards.h
Param.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Param.o: /usr/include/stdio.h /usr/include/limits.h /usr/include/sys/limits.h
Param.o: /usr/lpp/xlC/include/math.h /usr/include/math.h /usr/include/errno.h
Param.o: /usr/include/time.h types.h precision.h /usr/include/float.h
Param.o: Complex.h utils.h /usr/lpp/xlC/include/stddef.h
Param.o: /usr/include/stddef.h /usr/lpp/xlC/include/stdarg.h
Param.o: /usr/include/stdarg.h /usr/include/va_list.h
Param.o: /usr/lpp/xlC/include/stdlib.h /usr/include/stdlib.h new.h minmax.h
Param.o: pow.h dirsep.h DynVector.h Table.h Param.h Integrator.h Approx.h
ThreeWave.o: NWave.h kernel.h /usr/lpp/xlC/include/iostream.h
ThreeWave.o: /usr/include/memory.h /usr/include/string.h
ThreeWave.o: /usr/include/standards.h /usr/include/sys/types.h
ThreeWave.o: /usr/lpp/xlC/include/fstream.h /usr/include/stdio.h
ThreeWave.o: /usr/include/limits.h /usr/include/sys/limits.h
ThreeWave.o: /usr/lpp/xlC/include/math.h /usr/include/math.h
ThreeWave.o: /usr/include/errno.h /usr/include/time.h types.h precision.h
ThreeWave.o: /usr/include/float.h Complex.h utils.h
ThreeWave.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
ThreeWave.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
ThreeWave.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
ThreeWave.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
ThreeWave.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Navier.o: NWave.h kernel.h /usr/lpp/xlC/include/iostream.h
Navier.o: /usr/include/memory.h /usr/include/string.h
Navier.o: /usr/include/standards.h /usr/include/sys/types.h
Navier.o: /usr/lpp/xlC/include/fstream.h /usr/include/stdio.h
Navier.o: /usr/include/limits.h /usr/include/sys/limits.h
Navier.o: /usr/lpp/xlC/include/math.h /usr/include/math.h
Navier.o: /usr/include/errno.h /usr/include/time.h types.h precision.h
Navier.o: /usr/include/float.h Complex.h utils.h
Navier.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
Navier.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
Navier.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
Navier.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
Navier.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Navier.o: Polar.h
NWave.o: NWave.h kernel.h /usr/lpp/xlC/include/iostream.h
NWave.o: /usr/include/memory.h /usr/include/string.h /usr/include/standards.h
NWave.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
NWave.o: /usr/include/stdio.h /usr/include/limits.h /usr/include/sys/limits.h
NWave.o: /usr/lpp/xlC/include/math.h /usr/include/math.h /usr/include/errno.h
NWave.o: /usr/include/time.h types.h precision.h /usr/include/float.h
NWave.o: Complex.h utils.h /usr/lpp/xlC/include/stddef.h
NWave.o: /usr/include/stddef.h /usr/lpp/xlC/include/stdarg.h
NWave.o: /usr/include/stdarg.h /usr/include/va_list.h
NWave.o: /usr/lpp/xlC/include/stdlib.h /usr/include/stdlib.h new.h minmax.h
NWave.o: pow.h dirsep.h DynVector.h Table.h Param.h Integrator.h Approx.h
NWave.o: Geometry.h Pair.h Bin.h
Geometry.o: NWave.h kernel.h /usr/lpp/xlC/include/iostream.h
Geometry.o: /usr/include/memory.h /usr/include/string.h
Geometry.o: /usr/include/standards.h /usr/include/sys/types.h
Geometry.o: /usr/lpp/xlC/include/fstream.h /usr/include/stdio.h
Geometry.o: /usr/include/limits.h /usr/include/sys/limits.h
Geometry.o: /usr/lpp/xlC/include/math.h /usr/include/math.h
Geometry.o: /usr/include/errno.h /usr/include/time.h types.h precision.h
Geometry.o: /usr/include/float.h Complex.h utils.h
Geometry.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
Geometry.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
Geometry.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
Geometry.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
Geometry.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Polar.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Polar.o: /usr/include/string.h /usr/include/standards.h
Polar.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Polar.o: /usr/include/stdio.h /usr/include/limits.h /usr/include/sys/limits.h
Polar.o: /usr/lpp/xlC/include/math.h /usr/include/math.h /usr/include/errno.h
Polar.o: /usr/include/time.h types.h precision.h /usr/include/float.h
Polar.o: Complex.h utils.h /usr/lpp/xlC/include/stddef.h
Polar.o: /usr/include/stddef.h /usr/lpp/xlC/include/stdarg.h
Polar.o: /usr/include/stdarg.h /usr/include/va_list.h
Polar.o: /usr/lpp/xlC/include/stdlib.h /usr/include/stdlib.h new.h minmax.h
Polar.o: pow.h dirsep.h DynVector.h Table.h Param.h Integrator.h Approx.h
Polar.o: Geometry.h Pair.h Bin.h Polar.h
PolarAverage.o: Polar.h types.h precision.h /usr/include/float.h Complex.h
PolarAverage.o: /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
PolarAverage.o: /usr/include/string.h /usr/include/standards.h
PolarAverage.o: /usr/include/sys/types.h /usr/lpp/xlC/include/math.h
PolarAverage.o: /usr/include/math.h Bin.h utils.h
PolarAverage.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
PolarAverage.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
PolarAverage.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
PolarAverage.o: /usr/include/stdlib.h /usr/lpp/xlC/include/fstream.h new.h
PolarAverage.o: minmax.h pow.h /usr/include/limits.h
PolarAverage.o: /usr/include/sys/limits.h dirsep.h
simpfast.o: /usr/lpp/xlC/include/math.h /usr/include/math.h precision.h
simpfast.o: /usr/include/float.h
Kepler.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Kepler.o: /usr/include/string.h /usr/include/standards.h
Kepler.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Kepler.o: /usr/include/stdio.h /usr/include/limits.h
Kepler.o: /usr/include/sys/limits.h /usr/lpp/xlC/include/math.h
Kepler.o: /usr/include/math.h /usr/include/errno.h /usr/include/time.h
Kepler.o: types.h precision.h /usr/include/float.h Complex.h utils.h
Kepler.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
Kepler.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
Kepler.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
Kepler.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
Kepler.o: Table.h Param.h Integrator.h Approx.h
Lotka.o: kernel.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
Lotka.o: /usr/include/string.h /usr/include/standards.h
Lotka.o: /usr/include/sys/types.h /usr/lpp/xlC/include/fstream.h
Lotka.o: /usr/include/stdio.h /usr/include/limits.h /usr/include/sys/limits.h
Lotka.o: /usr/lpp/xlC/include/math.h /usr/include/math.h /usr/include/errno.h
Lotka.o: /usr/include/time.h types.h precision.h /usr/include/float.h
Lotka.o: Complex.h utils.h /usr/lpp/xlC/include/stddef.h
Lotka.o: /usr/include/stddef.h /usr/lpp/xlC/include/stdarg.h
Lotka.o: /usr/include/stdarg.h /usr/include/va_list.h
Lotka.o: /usr/lpp/xlC/include/stdlib.h /usr/include/stdlib.h new.h minmax.h
Lotka.o: pow.h dirsep.h DynVector.h Table.h Param.h Integrator.h Approx.h
utils.o: /usr/include/ctype.h /usr/include/standards.h
utils.o: /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
utils.o: /usr/include/string.h /usr/include/sys/types.h /usr/include/stdio.h
utils.o: /usr/lpp/xlC/include/unistd.h /usr/include/unistd.h
utils.o: /usr/include/sys/access.h /usr/include/errno.h
utils.o: /usr/include/sys/times.h /usr/include/time.h utils.h
utils.o: /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
utils.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
utils.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
utils.o: /usr/include/stdlib.h /usr/lpp/xlC/include/fstream.h
utils.o: /usr/lpp/xlC/include/math.h /usr/include/math.h new.h precision.h
utils.o: /usr/include/float.h minmax.h Complex.h pow.h /usr/include/limits.h
utils.o: /usr/include/sys/limits.h dirsep.h
strcasecmp.o: /usr/include/string.h /usr/include/standards.h
strcasecmp.o: /usr/include/sys/types.h
new.o: new.h /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
new.o: /usr/include/string.h /usr/include/standards.h
new.o: /usr/include/sys/types.h
idle.o: /usr/lpp/xlC/include/unistd.h /usr/include/unistd.h
idle.o: /usr/include/sys/access.h /usr/include/standards.h kernel.h
idle.o: /usr/lpp/xlC/include/iostream.h /usr/include/memory.h
idle.o: /usr/include/string.h /usr/include/sys/types.h
idle.o: /usr/lpp/xlC/include/fstream.h /usr/include/stdio.h
idle.o: /usr/include/limits.h /usr/include/sys/limits.h
idle.o: /usr/lpp/xlC/include/math.h /usr/include/math.h /usr/include/errno.h
idle.o: /usr/include/time.h types.h precision.h /usr/include/float.h
idle.o: Complex.h utils.h /usr/lpp/xlC/include/stddef.h /usr/include/stddef.h
idle.o: /usr/lpp/xlC/include/stdarg.h /usr/include/stdarg.h
idle.o: /usr/include/va_list.h /usr/lpp/xlC/include/stdlib.h
idle.o: /usr/include/stdlib.h new.h minmax.h pow.h dirsep.h DynVector.h
idle.o: Table.h Param.h Integrator.h Approx.h
