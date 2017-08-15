lib_LTLIBRARIES = libsciutils.la

libsciutils_la_CPPFLAGS = -I.
#libsciutils_la_DEPENDENCIES = Makefile

if CFITSIO

libsciutils_la_SOURCES = sciutils.h genbetai.c readdata.c readxray.c ztompc.c alloc.c
libsciutils_la_CFLAGS = $(OPTFLAGS)

else

libsciutils_la_SOURCES = sciutils.h genbetai.c readdata.c ztompc.c alloc.c
libsciutils_la_CFLAGS = $(OPTFLAGS)

endif

noinst_LIBRARIES = libsciutils.a

libsciutils_a_CPPFLAGS = -I.
#libsciutils_a_DEPENDENCIES = Makefile

if CFITSIO

libsciutils_a_SOURCES = sciutils.h genbetai.c readdata.c readxray.c ztompc.c alloc.c
libsciutils_a_CFLAGS = $(OPTFLAGS)

else

libsciutils_a_SOURCES = sciutils.h genbetai.c readdata.c ztompc.c alloc.c
libsciutils_a_CFLAGS = $(OPTFLAGS)

endif

bin_PROGRAMS = angdist

angdist_SOURCES = angdist.c

angdist_LDADD = libsciutils.a $(LDFLAGS)

info_TEXINFOS = sciutils.texinfo

man3_MANS = sciutils.3

BUILT_SOURCES = sciutils.3
EXTRA_DIST = buildscripts/texi2man sciutils.pdf 

sciutils.3: sciutils.texinfo 
	buildscripts/texi2man sciutils.texinfo > $@
