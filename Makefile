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

kernel.o: /usr/lib/g++-include/iomanip.h /usr/lib/g++-include/iostream.h
kernel.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
kernel.o: /usr/lib/g++-include/_G_config.h kernel.h
kernel.o: /usr/lib/g++-include/fstream.h /usr/include/stdio.h
kernel.o: /usr/local/lib/gcc-include/limits.h
kernel.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
kernel.o: /usr/include/features.h /usr/include/sys/cdefs.h
kernel.o: /usr/include/huge_val.h /usr/include/endian.h
kernel.o: /usr/include/bytesex.h /usr/include/linux/version.h
kernel.o: /usr/include/asm/byteorder.h /usr/include/nan.h
kernel.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
kernel.o: /usr/local/include/i386/__math.h
kernel.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
kernel.o: /usr/include/linux/errno.h /usr/include/time.h
kernel.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
kernel.o: /usr/lib/g++-include/std/cstring.h
kernel.o: /usr/lib/g++-include/std/cstddef.h
kernel.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
kernel.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
kernel.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
kernel.o: Table.h Param.h Integrator.h Approx.h
Approx.o: kernel.h /usr/lib/g++-include/iostream.h
Approx.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Approx.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Approx.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Approx.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Approx.o: /usr/include/features.h /usr/include/sys/cdefs.h
Approx.o: /usr/include/huge_val.h /usr/include/endian.h
Approx.o: /usr/include/bytesex.h /usr/include/linux/version.h
Approx.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Approx.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Approx.o: /usr/local/include/i386/__math.h
Approx.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Approx.o: /usr/include/linux/errno.h /usr/include/time.h
Approx.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Approx.o: /usr/lib/g++-include/std/cstring.h
Approx.o: /usr/lib/g++-include/std/cstddef.h
Approx.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
Approx.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Approx.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Approx.o: Table.h Param.h Integrator.h Approx.h
Integrator.o: kernel.h /usr/lib/g++-include/iostream.h
Integrator.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Integrator.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Integrator.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Integrator.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Integrator.o: /usr/include/features.h /usr/include/sys/cdefs.h
Integrator.o: /usr/include/huge_val.h /usr/include/endian.h
Integrator.o: /usr/include/bytesex.h /usr/include/linux/version.h
Integrator.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Integrator.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Integrator.o: /usr/local/include/i386/__math.h
Integrator.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Integrator.o: /usr/include/linux/errno.h /usr/include/time.h
Integrator.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Integrator.o: /usr/lib/g++-include/std/cstring.h
Integrator.o: /usr/lib/g++-include/std/cstddef.h
Integrator.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
Integrator.o: utils.h /usr/local/lib/gcc-include/stdarg.h
Integrator.o: /usr/include/stdlib.h /usr/include/alloca.h new.h minmax.h
Integrator.o: pow.h dirsep.h DynVector.h Table.h Param.h Integrator.h
Integrator.o: Approx.h
Param.o: kernel.h /usr/lib/g++-include/iostream.h
Param.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Param.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Param.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Param.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Param.o: /usr/include/features.h /usr/include/sys/cdefs.h
Param.o: /usr/include/huge_val.h /usr/include/endian.h /usr/include/bytesex.h
Param.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Param.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Param.o: /usr/include/values.h /usr/local/include/i386/__math.h
Param.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Param.o: /usr/include/linux/errno.h /usr/include/time.h
Param.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Param.o: /usr/lib/g++-include/std/cstring.h
Param.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Param.o: types.h precision.h Complex.h utils.h
Param.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Param.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Param.o: Table.h Param.h Integrator.h Approx.h
ThreeWave.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
ThreeWave.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
ThreeWave.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
ThreeWave.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
ThreeWave.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
ThreeWave.o: /usr/include/features.h /usr/include/sys/cdefs.h
ThreeWave.o: /usr/include/huge_val.h /usr/include/endian.h
ThreeWave.o: /usr/include/bytesex.h /usr/include/linux/version.h
ThreeWave.o: /usr/include/asm/byteorder.h /usr/include/nan.h
ThreeWave.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
ThreeWave.o: /usr/local/include/i386/__math.h
ThreeWave.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
ThreeWave.o: /usr/include/linux/errno.h /usr/include/time.h
ThreeWave.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
ThreeWave.o: /usr/lib/g++-include/std/cstring.h
ThreeWave.o: /usr/lib/g++-include/std/cstddef.h
ThreeWave.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
ThreeWave.o: utils.h /usr/local/lib/gcc-include/stdarg.h
ThreeWave.o: /usr/include/stdlib.h /usr/include/alloca.h new.h minmax.h pow.h
ThreeWave.o: dirsep.h DynVector.h Table.h Param.h Integrator.h Approx.h
ThreeWave.o: Geometry.h Pair.h Bin.h
Navier.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
Navier.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Navier.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Navier.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Navier.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Navier.o: /usr/include/features.h /usr/include/sys/cdefs.h
Navier.o: /usr/include/huge_val.h /usr/include/endian.h
Navier.o: /usr/include/bytesex.h /usr/include/linux/version.h
Navier.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Navier.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Navier.o: /usr/local/include/i386/__math.h
Navier.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Navier.o: /usr/include/linux/errno.h /usr/include/time.h
Navier.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Navier.o: /usr/lib/g++-include/std/cstring.h
Navier.o: /usr/lib/g++-include/std/cstddef.h
Navier.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
Navier.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Navier.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Navier.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Navier.o: Polar.h
NWave.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
NWave.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
NWave.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
NWave.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
NWave.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
NWave.o: /usr/include/features.h /usr/include/sys/cdefs.h
NWave.o: /usr/include/huge_val.h /usr/include/endian.h /usr/include/bytesex.h
NWave.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
NWave.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
NWave.o: /usr/include/values.h /usr/local/include/i386/__math.h
NWave.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
NWave.o: /usr/include/linux/errno.h /usr/include/time.h
NWave.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
NWave.o: /usr/lib/g++-include/std/cstring.h
NWave.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
NWave.o: types.h precision.h Complex.h utils.h
NWave.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
NWave.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
NWave.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Geometry.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
Geometry.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Geometry.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Geometry.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Geometry.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Geometry.o: /usr/include/features.h /usr/include/sys/cdefs.h
Geometry.o: /usr/include/huge_val.h /usr/include/endian.h
Geometry.o: /usr/include/bytesex.h /usr/include/linux/version.h
Geometry.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Geometry.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Geometry.o: /usr/local/include/i386/__math.h
Geometry.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Geometry.o: /usr/include/linux/errno.h /usr/include/time.h
Geometry.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Geometry.o: /usr/lib/g++-include/std/cstring.h
Geometry.o: /usr/lib/g++-include/std/cstddef.h
Geometry.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
Geometry.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Geometry.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Geometry.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Polar.o: kernel.h /usr/lib/g++-include/iostream.h
Polar.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Polar.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Polar.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Polar.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Polar.o: /usr/include/features.h /usr/include/sys/cdefs.h
Polar.o: /usr/include/huge_val.h /usr/include/endian.h /usr/include/bytesex.h
Polar.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Polar.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Polar.o: /usr/include/values.h /usr/local/include/i386/__math.h
Polar.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Polar.o: /usr/include/linux/errno.h /usr/include/time.h
Polar.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Polar.o: /usr/lib/g++-include/std/cstring.h
Polar.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Polar.o: types.h precision.h Complex.h utils.h
Polar.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Polar.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Polar.o: Table.h Param.h Integrator.h Approx.h Geometry.h Pair.h Bin.h
Polar.o: Polar.h
PolarAverage.o: Polar.h types.h precision.h
PolarAverage.o: /usr/local/lib/gcc-include/float.h Complex.h
PolarAverage.o: /usr/lib/g++-include/iostream.h
PolarAverage.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
PolarAverage.o: /usr/lib/g++-include/_G_config.h /usr/include/math.h
PolarAverage.o: /usr/include/features.h /usr/include/sys/cdefs.h
PolarAverage.o: /usr/include/huge_val.h /usr/include/endian.h
PolarAverage.o: /usr/include/bytesex.h /usr/include/linux/version.h
PolarAverage.o: /usr/include/asm/byteorder.h /usr/include/nan.h
PolarAverage.o: /usr/include/values.h /usr/local/include/i386/__math.h
PolarAverage.o: /usr/local/include/i386/extensions.h Bin.h utils.h
PolarAverage.o: /usr/lib/g++-include/std/stddef.h
PolarAverage.o: /usr/lib/g++-include/std/cstddef.h
PolarAverage.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
PolarAverage.o: /usr/include/errno.h /usr/include/linux/errno.h
PolarAverage.o: /usr/include/alloca.h /usr/lib/g++-include/fstream.h new.h
PolarAverage.o: minmax.h pow.h /usr/local/lib/gcc-include/limits.h
PolarAverage.o: /usr/local/lib/gcc-include/syslimits.h dirsep.h
simpfast.o: /usr/include/math.h /usr/include/features.h
simpfast.o: /usr/include/sys/cdefs.h /usr/include/huge_val.h
simpfast.o: /usr/include/endian.h /usr/include/bytesex.h
simpfast.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
simpfast.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
simpfast.o: /usr/include/values.h /usr/local/include/i386/__math.h
simpfast.o: /usr/local/include/i386/extensions.h precision.h
Kepler.o: kernel.h /usr/lib/g++-include/iostream.h
Kepler.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Kepler.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Kepler.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Kepler.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Kepler.o: /usr/include/features.h /usr/include/sys/cdefs.h
Kepler.o: /usr/include/huge_val.h /usr/include/endian.h
Kepler.o: /usr/include/bytesex.h /usr/include/linux/version.h
Kepler.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Kepler.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Kepler.o: /usr/local/include/i386/__math.h
Kepler.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Kepler.o: /usr/include/linux/errno.h /usr/include/time.h
Kepler.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Kepler.o: /usr/lib/g++-include/std/cstring.h
Kepler.o: /usr/lib/g++-include/std/cstddef.h
Kepler.o: /usr/lib/g++-include/std/stddef.h types.h precision.h Complex.h
Kepler.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Kepler.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Kepler.o: Table.h Param.h Integrator.h Approx.h
Lotka.o: kernel.h /usr/lib/g++-include/iostream.h
Lotka.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Lotka.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Lotka.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Lotka.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
Lotka.o: /usr/include/features.h /usr/include/sys/cdefs.h
Lotka.o: /usr/include/huge_val.h /usr/include/endian.h /usr/include/bytesex.h
Lotka.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Lotka.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Lotka.o: /usr/include/values.h /usr/local/include/i386/__math.h
Lotka.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
Lotka.o: /usr/include/linux/errno.h /usr/include/time.h
Lotka.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Lotka.o: /usr/lib/g++-include/std/cstring.h
Lotka.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Lotka.o: types.h precision.h Complex.h utils.h
Lotka.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Lotka.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
Lotka.o: Table.h Param.h Integrator.h Approx.h
utils.o: /usr/include/ctype.h /usr/include/features.h
utils.o: /usr/include/sys/cdefs.h /usr/include/endian.h
utils.o: /usr/include/bytesex.h /usr/include/linux/version.h
utils.o: /usr/include/asm/byteorder.h /usr/lib/g++-include/iostream.h
utils.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
utils.o: /usr/lib/g++-include/_G_config.h /usr/include/stdio.h
utils.o: /usr/include/unistd.h /usr/include/posix_opt.h
utils.o: /usr/include/gnu/types.h /usr/lib/g++-include/std/stddef.h
utils.o: /usr/lib/g++-include/std/cstddef.h /usr/include/confname.h
utils.o: /usr/include/sys/types.h /usr/include/linux/types.h
utils.o: /usr/include/asm/types.h /usr/include/errno.h
utils.o: /usr/include/linux/errno.h /usr/include/sys/times.h
utils.o: /usr/include/time.h /usr/include/linux/times.h
utils.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
utils.o: /usr/lib/g++-include/std/cstring.h utils.h
utils.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
utils.o: /usr/include/alloca.h /usr/lib/g++-include/fstream.h
utils.o: /usr/include/math.h /usr/include/huge_val.h /usr/include/nan.h
utils.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
utils.o: /usr/local/include/i386/__math.h
utils.o: /usr/local/include/i386/extensions.h new.h precision.h minmax.h
utils.o: Complex.h pow.h /usr/local/lib/gcc-include/limits.h
utils.o: /usr/local/lib/gcc-include/syslimits.h dirsep.h
strcasecmp.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
strcasecmp.o: /usr/lib/g++-include/std/cstring.h
strcasecmp.o: /usr/lib/g++-include/std/cstddef.h
strcasecmp.o: /usr/lib/g++-include/std/stddef.h
strcasecmp.o: /usr/lib/g++-include/_G_config.h
new.o: /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h
new.o: /usr/lib/g++-include/std/stddef.h /usr/lib/g++-include/_G_config.h
new.o: /usr/lib/g++-include/std/cstddef.h /usr/include/errno.h
new.o: /usr/include/linux/errno.h /usr/include/alloca.h
new.o: /usr/lib/g++-include/iostream.h /usr/lib/g++-include/streambuf.h
new.o: /usr/lib/g++-include/libio.h new.h
idle.o: /usr/include/unistd.h /usr/include/features.h
idle.o: /usr/include/sys/cdefs.h /usr/include/posix_opt.h
idle.o: /usr/include/gnu/types.h /usr/lib/g++-include/std/stddef.h
idle.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/std/cstddef.h
idle.o: /usr/include/confname.h /usr/include/sys/types.h
idle.o: /usr/include/linux/types.h /usr/include/asm/types.h kernel.h
idle.o: /usr/lib/g++-include/iostream.h /usr/lib/g++-include/streambuf.h
idle.o: /usr/lib/g++-include/libio.h /usr/lib/g++-include/fstream.h
idle.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
idle.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/math.h
idle.o: /usr/include/huge_val.h /usr/include/endian.h /usr/include/bytesex.h
idle.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
idle.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
idle.o: /usr/include/values.h /usr/local/include/i386/__math.h
idle.o: /usr/local/include/i386/extensions.h /usr/include/errno.h
idle.o: /usr/include/linux/errno.h /usr/include/time.h
idle.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
idle.o: /usr/lib/g++-include/std/cstring.h types.h precision.h Complex.h
idle.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
idle.o: /usr/include/alloca.h new.h minmax.h pow.h dirsep.h DynVector.h
idle.o: Table.h Param.h Integrator.h Approx.h
