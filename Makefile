LIB = $(LFLAGS) -lm
OPT = -g $(CFLAGS)
MAKEDEPEND = makedepend

ARCH = unix.o
POLL = poll.o

include config/$(HOSTTYPE)

UTILS = utils.o strcasecmp.o new.o $(POLL) $(ARCH)
CORE = kernel.o Approx.o Integrator.o Param.o $(UTILS)
POLAR = Polar.o PolarAverage.o simpfast.o
TRIAD = $(CORE) NWave.o Geometry.o $(POLAR)


.SUFFIXES: .cc

triad:	Navier.o $(TRIAD)
	$(C++) $(OPT) -o triad Navier.o $(TRIAD) $(LIB) 

w3:	ThreeWave.o $(TRIAD)
	$(C++) $(OPT) -o triad ThreeWave.o NWave.o $(CORE) $(LIB)

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
	$(MAKEDEPEND) $(MDOPT) -I /usr/local/include \
	kernel.cc Approx.cc Integrator.cc Param.cc ThreeWave.cc \
	Navier.cc NWave.cc Geometry.cc Polar.cc PolarAverage.cc simpfast.cc \
	Kepler.cc Lotka.cc utils.cc strcasecmp.cc new.cc poll.cc \
	idle.cc unix.cc

# DO NOT DELETE THIS LINE -- make depend depends on it.

kernel.o: kernel.h /usr/lib/g++-include/iostream.h
kernel.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
kernel.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
kernel.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
kernel.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
kernel.o: /usr/include/features.h /usr/include/sys/cdefs.h
kernel.o: /usr/include/linux/errno.h /usr/include/time.h
kernel.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
kernel.o: /usr/lib/g++-include/std/cstring.h
kernel.o: /usr/lib/g++-include/std/cstddef.h
kernel.o: /usr/lib/g++-include/std/stddef.h utils.h
kernel.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
kernel.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
kernel.o: /usr/include/endian.h /usr/include/bytesex.h
kernel.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
kernel.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
kernel.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
kernel.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
kernel.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
kernel.o: /usr/lib/g++-include/iomanip.h
Approx.o: kernel.h /usr/lib/g++-include/iostream.h
Approx.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Approx.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Approx.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Approx.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Approx.o: /usr/include/features.h /usr/include/sys/cdefs.h
Approx.o: /usr/include/linux/errno.h /usr/include/time.h
Approx.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Approx.o: /usr/lib/g++-include/std/cstring.h
Approx.o: /usr/lib/g++-include/std/cstddef.h
Approx.o: /usr/lib/g++-include/std/stddef.h utils.h
Approx.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Approx.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Approx.o: /usr/include/endian.h /usr/include/bytesex.h
Approx.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Approx.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Approx.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Approx.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Approx.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
Integrator.o: kernel.h /usr/lib/g++-include/iostream.h
Integrator.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Integrator.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Integrator.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Integrator.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Integrator.o: /usr/include/features.h /usr/include/sys/cdefs.h
Integrator.o: /usr/include/linux/errno.h /usr/include/time.h
Integrator.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Integrator.o: /usr/lib/g++-include/std/cstring.h
Integrator.o: /usr/lib/g++-include/std/cstddef.h
Integrator.o: /usr/lib/g++-include/std/stddef.h utils.h
Integrator.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Integrator.o: /usr/include/alloca.h /usr/include/math.h
Integrator.o: /usr/include/huge_val.h /usr/include/endian.h
Integrator.o: /usr/include/bytesex.h /usr/include/linux/version.h
Integrator.o: /usr/include/asm/byteorder.h /usr/include/nan.h
Integrator.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
Integrator.o: /usr/local/include/i386/__math.h arch/i386.h
Integrator.o: /usr/local/include/i386/extensions.h new.h precision.h
Integrator.o: Complex.h pow.h types.h DynVector.h Table.h Param.h
Integrator.o: Integrator.h Approx.h
Param.o: kernel.h /usr/lib/g++-include/iostream.h
Param.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Param.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Param.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Param.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Param.o: /usr/include/features.h /usr/include/sys/cdefs.h
Param.o: /usr/include/linux/errno.h /usr/include/time.h
Param.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Param.o: /usr/lib/g++-include/std/cstring.h
Param.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Param.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Param.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Param.o: /usr/include/endian.h /usr/include/bytesex.h
Param.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Param.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Param.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Param.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Param.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
ThreeWave.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
ThreeWave.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
ThreeWave.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
ThreeWave.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
ThreeWave.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
ThreeWave.o: /usr/include/features.h /usr/include/sys/cdefs.h
ThreeWave.o: /usr/include/linux/errno.h /usr/include/time.h
ThreeWave.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
ThreeWave.o: /usr/lib/g++-include/std/cstring.h
ThreeWave.o: /usr/lib/g++-include/std/cstddef.h
ThreeWave.o: /usr/lib/g++-include/std/stddef.h utils.h
ThreeWave.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
ThreeWave.o: /usr/include/alloca.h /usr/include/math.h
ThreeWave.o: /usr/include/huge_val.h /usr/include/endian.h
ThreeWave.o: /usr/include/bytesex.h /usr/include/linux/version.h
ThreeWave.o: /usr/include/asm/byteorder.h /usr/include/nan.h
ThreeWave.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
ThreeWave.o: /usr/local/include/i386/__math.h arch/i386.h
ThreeWave.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
ThreeWave.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
ThreeWave.o: Geometry.h Pair.h Bin.h
Navier.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
Navier.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Navier.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Navier.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Navier.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Navier.o: /usr/include/features.h /usr/include/sys/cdefs.h
Navier.o: /usr/include/linux/errno.h /usr/include/time.h
Navier.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Navier.o: /usr/lib/g++-include/std/cstring.h
Navier.o: /usr/lib/g++-include/std/cstddef.h
Navier.o: /usr/lib/g++-include/std/stddef.h utils.h
Navier.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Navier.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Navier.o: /usr/include/endian.h /usr/include/bytesex.h
Navier.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Navier.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Navier.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Navier.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Navier.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
Navier.o: Geometry.h Pair.h Bin.h Polar.h
NWave.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
NWave.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
NWave.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
NWave.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
NWave.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
NWave.o: /usr/include/features.h /usr/include/sys/cdefs.h
NWave.o: /usr/include/linux/errno.h /usr/include/time.h
NWave.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
NWave.o: /usr/lib/g++-include/std/cstring.h
NWave.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
NWave.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
NWave.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
NWave.o: /usr/include/endian.h /usr/include/bytesex.h
NWave.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
NWave.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
NWave.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
NWave.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
NWave.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
NWave.o: Geometry.h Pair.h Bin.h
Geometry.o: NWave.h kernel.h /usr/lib/g++-include/iostream.h
Geometry.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Geometry.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Geometry.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Geometry.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Geometry.o: /usr/include/features.h /usr/include/sys/cdefs.h
Geometry.o: /usr/include/linux/errno.h /usr/include/time.h
Geometry.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Geometry.o: /usr/lib/g++-include/std/cstring.h
Geometry.o: /usr/lib/g++-include/std/cstddef.h
Geometry.o: /usr/lib/g++-include/std/stddef.h utils.h
Geometry.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Geometry.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Geometry.o: /usr/include/endian.h /usr/include/bytesex.h
Geometry.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Geometry.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Geometry.o: /usr/include/values.h /usr/local/include/i386/__math.h
Geometry.o: arch/i386.h /usr/local/include/i386/extensions.h new.h
Geometry.o: precision.h Complex.h pow.h types.h DynVector.h Table.h Param.h
Geometry.o: Integrator.h Approx.h Geometry.h Pair.h Bin.h
Polar.o: kernel.h /usr/lib/g++-include/iostream.h
Polar.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Polar.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Polar.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Polar.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Polar.o: /usr/include/features.h /usr/include/sys/cdefs.h
Polar.o: /usr/include/linux/errno.h /usr/include/time.h
Polar.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Polar.o: /usr/lib/g++-include/std/cstring.h
Polar.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Polar.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Polar.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Polar.o: /usr/include/endian.h /usr/include/bytesex.h
Polar.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Polar.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Polar.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Polar.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Polar.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
Polar.o: Geometry.h Pair.h Bin.h Polar.h
PolarAverage.o: Polar.h types.h precision.h
PolarAverage.o: /usr/local/lib/gcc-include/float.h Complex.h
PolarAverage.o: /usr/lib/g++-include/iostream.h
PolarAverage.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
PolarAverage.o: /usr/lib/g++-include/_G_config.h /usr/include/math.h
PolarAverage.o: /usr/include/features.h /usr/include/sys/cdefs.h
PolarAverage.o: /usr/include/huge_val.h /usr/include/endian.h
PolarAverage.o: /usr/include/bytesex.h /usr/include/linux/version.h
PolarAverage.o: /usr/include/asm/byteorder.h /usr/include/nan.h
PolarAverage.o: /usr/include/values.h /usr/local/include/i386/__math.h Bin.h
PolarAverage.o: utils.h /usr/lib/g++-include/std/stddef.h
PolarAverage.o: /usr/lib/g++-include/std/cstddef.h
PolarAverage.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
PolarAverage.o: /usr/include/errno.h /usr/include/linux/errno.h
PolarAverage.o: /usr/include/alloca.h /usr/lib/g++-include/fstream.h
PolarAverage.o: arch/i386.h /usr/local/include/i386/extensions.h new.h pow.h
PolarAverage.o: /usr/local/lib/gcc-include/limits.h
PolarAverage.o: /usr/local/lib/gcc-include/syslimits.h Geometry.h kernel.h
PolarAverage.o: /usr/include/stdio.h /usr/include/time.h
PolarAverage.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
PolarAverage.o: /usr/lib/g++-include/std/cstring.h DynVector.h Table.h
PolarAverage.o: Param.h Integrator.h Approx.h Pair.h
simpfast.o: /usr/include/math.h /usr/include/features.h
simpfast.o: /usr/include/sys/cdefs.h /usr/include/huge_val.h
simpfast.o: /usr/include/endian.h /usr/include/bytesex.h
simpfast.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
simpfast.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
simpfast.o: /usr/include/values.h /usr/local/include/i386/__math.h
simpfast.o: precision.h
Kepler.o: kernel.h /usr/lib/g++-include/iostream.h
Kepler.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Kepler.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Kepler.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Kepler.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Kepler.o: /usr/include/features.h /usr/include/sys/cdefs.h
Kepler.o: /usr/include/linux/errno.h /usr/include/time.h
Kepler.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Kepler.o: /usr/lib/g++-include/std/cstring.h
Kepler.o: /usr/lib/g++-include/std/cstddef.h
Kepler.o: /usr/lib/g++-include/std/stddef.h utils.h
Kepler.o: /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Kepler.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Kepler.o: /usr/include/endian.h /usr/include/bytesex.h
Kepler.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Kepler.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Kepler.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Kepler.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Kepler.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
Lotka.o: kernel.h /usr/lib/g++-include/iostream.h
Lotka.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
Lotka.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/fstream.h
Lotka.o: /usr/include/stdio.h /usr/local/lib/gcc-include/limits.h
Lotka.o: /usr/local/lib/gcc-include/syslimits.h /usr/include/errno.h
Lotka.o: /usr/include/features.h /usr/include/sys/cdefs.h
Lotka.o: /usr/include/linux/errno.h /usr/include/time.h
Lotka.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
Lotka.o: /usr/lib/g++-include/std/cstring.h
Lotka.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
Lotka.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
Lotka.o: /usr/include/alloca.h /usr/include/math.h /usr/include/huge_val.h
Lotka.o: /usr/include/endian.h /usr/include/bytesex.h
Lotka.o: /usr/include/linux/version.h /usr/include/asm/byteorder.h
Lotka.o: /usr/include/nan.h /usr/local/lib/gcc-include/float.h
Lotka.o: /usr/include/values.h /usr/local/include/i386/__math.h arch/i386.h
Lotka.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
Lotka.o: pow.h types.h DynVector.h Table.h Param.h Integrator.h Approx.h
utils.o: /usr/include/ctype.h /usr/include/features.h
utils.o: /usr/include/sys/cdefs.h /usr/include/endian.h
utils.o: /usr/include/bytesex.h /usr/include/linux/version.h
utils.o: /usr/include/asm/byteorder.h /usr/lib/g++-include/iostream.h
utils.o: /usr/lib/g++-include/streambuf.h /usr/lib/g++-include/libio.h
utils.o: /usr/lib/g++-include/_G_config.h /usr/include/stdio.h
utils.o: /usr/include/errno.h /usr/include/linux/errno.h
utils.o: /usr/lib/g++-include/string.h /usr/lib/g++-include/cstring
utils.o: /usr/lib/g++-include/std/cstring.h
utils.o: /usr/lib/g++-include/std/cstddef.h /usr/lib/g++-include/std/stddef.h
utils.o: utils.h /usr/local/lib/gcc-include/stdarg.h /usr/include/stdlib.h
utils.o: /usr/include/alloca.h /usr/lib/g++-include/fstream.h
utils.o: /usr/include/math.h /usr/include/huge_val.h /usr/include/nan.h
utils.o: /usr/local/lib/gcc-include/float.h /usr/include/values.h
utils.o: /usr/local/include/i386/__math.h arch/i386.h
utils.o: /usr/local/include/i386/extensions.h new.h precision.h Complex.h
utils.o: pow.h /usr/local/lib/gcc-include/limits.h
utils.o: /usr/local/lib/gcc-include/syslimits.h
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
idle.o: /usr/include/stdio.h /usr/lib/g++-include/libio.h
idle.o: /usr/lib/g++-include/_G_config.h /usr/include/stdlib.h
idle.o: /usr/include/features.h /usr/include/sys/cdefs.h
idle.o: /usr/lib/g++-include/std/stddef.h /usr/lib/g++-include/std/cstddef.h
idle.o: /usr/include/errno.h /usr/include/linux/errno.h /usr/include/alloca.h
idle.o: /usr/include/unistd.h /usr/include/posix_opt.h
idle.o: /usr/include/gnu/types.h /usr/include/confname.h
idle.o: /usr/include/sys/types.h /usr/include/linux/types.h
idle.o: /usr/include/asm/types.h /usr/lib/g++-include/string.h
idle.o: /usr/lib/g++-include/cstring /usr/lib/g++-include/std/cstring.h
unix.o: /usr/include/stdlib.h /usr/include/features.h
unix.o: /usr/include/sys/cdefs.h /usr/lib/g++-include/std/stddef.h
unix.o: /usr/lib/g++-include/_G_config.h /usr/lib/g++-include/std/cstddef.h
unix.o: /usr/include/errno.h /usr/include/linux/errno.h /usr/include/alloca.h
unix.o: /usr/include/unistd.h /usr/include/posix_opt.h
unix.o: /usr/include/gnu/types.h /usr/include/confname.h
unix.o: /usr/include/sys/types.h /usr/include/linux/types.h
unix.o: /usr/include/asm/types.h /usr/include/stdio.h
unix.o: /usr/lib/g++-include/libio.h /usr/include/pwd.h
unix.o: /usr/include/sys/times.h /usr/include/time.h
unix.o: /usr/include/linux/times.h /usr/lib/g++-include/string.h
unix.o: /usr/lib/g++-include/cstring /usr/lib/g++-include/std/cstring.h
