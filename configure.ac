#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(JACO, 1.3, amahdavi@sfsu.edu)
AC_CONFIG_SRCDIR([jaco.c])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADER([config.h])
AC_PROG_LIBTOOL

LDFLAGS="$LDFLAGS -L/opt/local/lib -L/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/ -L/usr/local/lib -L$HOME/lib"
OPTFLAGS=-O3
CFLAGS="-Wall"

# Uses acx_mpi.m4

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([gsl/gsl_version.h stdlib.h string.h strings.h unistd.h])

# Check for basic math library
AC_CHECK_LIB([m],[sqrt],
  ,
  AC_MSG_ERROR("Could not find basic C math library!"))

# Check for OpenMP
# Check for OpenMP flags
AC_CHECK_LIB([gomp],[omp_get_num_threads],OPENMPL="-lgomp -fopenmp",OPENMPL="")
AC_CHECK_LIB([gomp],[omp_get_num_threads],OPTFLAGS="$OPTFLAGS -fopenmp",
      AC_MSG_NOTICE([Not compiling OpenMP support (required for visualization)]))

LDFLAGS="$LDFLAGS $OPENMPL"

# Check for tree vectorization
AX_CHECK_COMPILER_FLAGS([-ftree-vectorize],
          OPTFLAGS="$OPTFLAGS -ftree-vectorize",	
          AC_MSG_NOTICE([Tree-vectorizing not supported by compiler]))

# Check for msse
case $host in
   powerpc*-*-*) OPTFLAGS="$OPTFLAGS -maltivec" ;;
   *) OPTFLAGS="$OPTFLAGS -mfpmath=sse -msse2" ;;
esac
AC_SUBST(OPTFLAGS,"$OPTFLAGS")

# Optional libraries
AC_ARG_WITH([mpi], 
  [AS_HELP_STRING([--without-mpi],
      [Disable MPI multiprocessing support])],
  [],
  [with_mpi=yes])

AS_IF([test "x$with_mpi" != xno],
   [ACX_MPI(,
      AC_MSG_ERROR([Couldn't find MPI libraries (use --without-mpi to disable)]))],
      AC_MSG_NOTICE([Not compiling MPI support (required for clustering).]))
AM_CONDITIONAL(MPI,test x$with_mpi = xyes)

AS_IF([test "x$with_mpi" != xno],,MPICC="gcc")

AC_ARG_WITH([cfitsio], 
  [AS_HELP_STRING([--without-cfitsio],
      [Disable cfitsio (and standalone JACO)])],
  [],
  [with_cfitsio=yes])

# Other miscellaneous options

AC_ARG_WITH([debug], 
  [AS_HELP_STRING([--with-debug],
      [Disable optimization/enable debug symbols])],
  [AC_SUBST(OPTFLAGS,"-O0 -g")],
  [with_debug=no])

AC_ARG_WITH([sz], 
  [AS_HELP_STRING([--with-sz],
      [Enable SZ support])],
  [],
  [with_sz=no])
AM_CONDITIONAL(SZ,test x$with_sz = xyes)

AS_IF([test "x$with_sz" != xno],
  [AC_DEFINE([WITH_SZ],1,[ Define 1 if compiling SZ support ])
   AC_CHECK_LIB([fftw3],[fftw_plan_dft_2d])])

AS_IF([test "x$with_sz" != xno],
   AC_CHECK_LIB([pthread],[sem_wait],,
	AC_MSG_ERROR([Couldn't fin threaded lpthread])))

AS_IF([test "x$with_sz" != xno],
   AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],,
	AC_MSG_ERROR([Couldn't find threaded libfftw3])))


AC_ARG_WITH([bolocam], 
  [AS_HELP_STRING([--with-bolocam],
      [Enable Bolocam support])],
  [],
  [with_bolocam=no])
AM_CONDITIONAL(BOLOCAM,test x$with_bolocam = xyes)

AS_IF([test "x$with_bolocam" != xno],
  [AC_DEFINE([WITH_BOLOCAM],1,[ Define 1 if compiling Bolocam support ])
   AC_CHECK_LIB([fftw3],[fftw_plan_dft_2d])])

AS_IF([test "x$with_bolocam" != xno],
   AC_CHECK_LIB([pthread],[sem_wait],,
	AC_MSG_ERROR([Couldn't fin threaded lpthread])))

AS_IF([test "x$with_bolocam" != xno],
   AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],,
	AC_MSG_ERROR([Couldn't find threaded libfftw3])))

AS_IF([test "x$with_sz" != xno],
   AS_IF([test "x$with_bolocam" != xno],
   	       AC_MSG_ERROR([Currently can't configure both Bolocam and interferometer support]),),)

AC_ARG_WITH([libpath], 
  [AS_HELP_STRING([--with-libpath],
      [Hardcode the install library path])],
  [],
  [with_libpath=no])
AM_CONDITIONAL(LIBPATH,test x$with_libpath = xyes)

LIBCFITSIO= 
AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_LIB([cfitsio],[ffgky],
  [AC_DEFINE([HAVE_CFITSIO],1,[ Define 1 if compiling CFITSIO support ])
   AC_SUBST(LIBCFITSIO,"-lcfitsio")],
   [AC_MSG_ERROR([Couldn't find cfitsio libraries. Without them, you will
	         only be able to run JACO under the CXC Sherpa package.
                 No standalone support or scripted data reduction 
                support is possible. It is recommended, but not
                required, that you install cfitsio.
                To compile without cfitsio use --without-cfitsio.])])],
  [AC_MSG_WARN("Not building cfitsio---no standalone functionality.")])
AM_CONDITIONAL(CFITSIO,test x$with_cfitsio = xyes)

AC_CONFIG_SUBDIRS([hrothgar])

AC_CHECK_LIB([gslcblas], [cblas_dnrm2],
  ,
  [AC_MSG_ERROR("Could not find GSL library gslcblas.")])

AC_CHECK_LIB([z],[gzwrite],
  ,
  AC_MSG_ERROR("Could not find libz compression library!"))

AC_CHECK_LIB([gsl], [gsl_vector_alloc],
  ,
  [AC_MSG_ERROR("Could not find the main GSL library.")])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow strstr])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
