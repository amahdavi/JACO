if CFITSIO
FITSCONF = hrothgar_image.c
else
FITSCONF =
endif

if MPI
MPIROUTINE = hrothgar_mpi.c
else
MPIROUTINE =
endif

AM_CFLAGS = -Wall -Wfloat-equal -W -Wshadow -Wbad-function-cast -Wredundant-decls -Wdisabled-optimization  -std=c99 -Wno-unused-parameter -Wno-unknown-pragmas $(OPTFLAGS)
AM_CPPFLAGS = -Ihrothgar -Ihrothgar/sciutils -I/usr/local/include -I/usr/local/include/openmpi -I/opt/local/include -I/opt/local/include/openmpi -I/opt/open-mpi/tcp-gnu42/include/ -I/usr/include/openmpi -I/usr/local/include -I$(HOME)/include -Isciutils
AM_LDFLAGS = -L/opt/local/lib -L/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/ -L/usr/local/lib -L$(HOME)/lib

SUBDIRS = sciutils

BUILT_SOURCES = hrothgar_help.h hrothgar.3 
EXTRA_DIST = buildscripts/texi2c buildscripts/texi2man hrothgar.pdf hrothgar-example.pdf

lib_LTLIBRARIES = libhrothgar.la

libhrothgar_la_SOURCES = hrothgar.h hrothgar_help.h hrothgar_proto.h hrothgar.c $(MPIROUTINE) $(FITSCONF)
libhrothgar_la_CFLAGS = $(AM_CFLAGS)

noinst_LIBRARIES = libhrothgar.a

libhrothgar_a_SOURCES = hrothgar.h hrothgar_help.h hrothgar_proto.h hrothgar.c $(MPIROUTINE) $(FITSCONF)
libhrothgar_a_LIBADD = 
libhrothgar_a_CFLAGS = $(AM_CFLAGS)

bin_PROGRAMS = hrothgar_test mcmcprob mcfit

hrothgar_test_SOURCES = hrothgar_test.c
hrothgar_test_LDADD = libhrothgar.a sciutils/libsciutils.a

hrothgar_test_LINK = $(MPICC) -o $@  $(LDFLAGS)
hrothgar_test_CFLAGS = $(AM_CFLAGS)

mcmcprob_SOURCES = mcmcprob.c
mcmcprob_LDADD = libhrothgar.a sciutils/libsciutils.a
mcmcprob_CXXFLAGS = $(AM_CFLAGS)

mcfit_SOURCES = mcfit.c 
mcfit_LDADD = libhrothgar.a sciutils/libsciutils.a 
mcfit_CFLAGS = $(AM_CFLAGS)
mcfit_LINK = $(MPICC) -o $@  $(LDFLAGS) 


man3_MANS = hrothgar.3

info_TEXINFOS = hrothgar.texinfo 
hrothgar_TEXINFOS = commandline.texinfo hrothgar_test.c

hrothgar_help.h: commandline.texinfo
	buildscripts/texi2c commandline.texinfo > $@

hrothgar.3: hrothgar.texinfo commandline.texinfo hrothgar_test.c
	buildscripts/texi2man hrothgar.texinfo > $@

