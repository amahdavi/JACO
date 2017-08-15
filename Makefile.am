if SZ
SZSOURCES = sz_fitting.c sz_fitting.h 
SZFLAGS = -Isz 
SZLIBS =  -Lsz -llapack -lcblas -lf77blas -latlas  -ltweaked_recipes -lrecipes -lf2c
else 
SZSOURCES = 
SZEXTRASOURCES = 
SZFLAGS =
SZLIBS = 
endif


if BOLOCAM
SZSOURCES += sz_jaco_interface.c sz_jaco_interface.h 
endif

AM_CFLAGS = -Wall -Wfloat-equal -W -Wshadow -Wbad-function-cast -Wredundant-decls -Wdisabled-optimization -std=gnu99 -Wno-unused-parameter -Wno-unknown-pragmas -funroll-loops $(OPTFLAGS)
AM_CPPFLAGS = -Ihrothgar -Ihrothgar/sciutils -I/usr/local/include -I/usr/local/include/openmpi -I/opt/local/include -I/opt/local/include/openmpi -I/opt/open-mpi/tcp-gnu42/include/ -I/usr/include/openmpi -I/usr/local/include -I$(HOME)/include -Isciutils
AM_LDFLAGS = -L/opt/local/lib -L/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/ -L/usr/local/lib -L$(HOME)/lib

if LIBPATH
ACCELERATE = -Wl,--rpath -Wl,$(libdir)
else
ACCELERATE = 
endif

SUBDIRS = hrothgar

if CFITSIO

bin_PROGRAMS = jaco jacong psf radial profile
jaco_SOURCES = jaco.c main.c jaco.h veldisp.c sherpa.c params.h standalone.c observesim.c standalone.h $(SZSOURCES)
jaco_LDADD =  $(ZLIB) hrothgar/libhrothgar.a hrothgar/sciutils/libsciutils.a $(MPILIBS) $(LIBCFITSIO) $(SZLIBS) 
jaco_LINK = $(MPICC) -o $@  $(LDFLAGS) 
jaco_CFLAGS = $(SZFLAGS) $(AM_CFLAGS)
jaco_LDFLAGS = -L$(libdir) $(ACCELERATE) 

jacong_SOURCES = jaco.c main.c jaco.h veldisp.c sherpa.c params.h standalone.c observesim.c standalone.h $(SZSOURCES)
jacong_LDADD =  $(ZLIB) hrothgar/libhrothgar.a hrothgar/sciutils/libsciutils.a $(MPILIBS) $(LIBCFITSIO) $(SZLIBS) 
jacong_CPPFLAGS = -DNONGRAVITY $(AM_CPPFLAGS) 
jacong_LINK = $(MPICC) -o $@  $(LDFLAGS) 
jacong_CFLAGS =  $(SZFLAGS) $(AM_CFLAGS)
jacong_LDFLAGS = -L$(libdir) $(ACCELERATE) 

radial_SOURCES = radial.c 
radial_LDADD = hrothgar/sciutils/libsciutils.a $(LIBCFITSIO)
radial_CFLAGS = $(AM_CFLAGS)
radial_LDFLAGS = -L$(libdir) $(ACCELERATE)

profile_SOURCES = profile.c 
profile_LDADD = $(ZLIB) hrothgar/libhrothgar.a hrothgar/sciutils/libsciutils.a $(LIBCFITSIO)
profile_CFLAGS = $(AM_CFLAGS)
profile_LDFLAGS = -L$(libdir) $(ACCELERATE)
profile_LINK = $(MPICC) -o $@  $(LDFLAGS) 

else

bin_PROGRAMS = psf

endif

psf_SOURCES = psf.c 
psf_LDADD = hrothgar/sciutils/libsciutils.a
psf_CFLAGS = $(AM_CFLAGS)
psf_LDFLAGS = -L$(libdir) $(ACCELERATE)

smdir=$(datadir)/sm
sm_DATA = qplot.sm specplot.sm

bin_SCRIPTS = mixprofmcmc specplot epstopdf fastng fastjaco
EXTRA_DIST = qplot.sm specplot.sm mixprofmcmc specplot epstopdf fastng fastjaco
