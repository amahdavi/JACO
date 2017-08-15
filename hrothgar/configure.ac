#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Hrothgar, 2.3, andisheh.mahdavi@gmail.com)
AC_CONFIG_SRCDIR([hrothgar.c])
AM_INIT_AUTOMAKE([gnu])
AC_CONFIG_HEADER([config.h])
AM_PROG_CC_C_O
AC_PROG_LIBTOOL

LDFLAGS="$LDFLAGS -L/opt/local/lib -L/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/ -L/usr/local/lib -L$HOME/lib"
CPPFLAGS="$CPPFLAGS -I/opt/local/include -I/opt/local/include/openmpi"
OPTFLAGS=-O3

# Uses acx_mpi.m4

# Checks for programs.
# AC_PROG_CC
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([gsl/gsl_version.h stdlib.h string.h strings.h unistd.h])

# Basic math library
AC_CHECK_LIB([m],[sqrt],
  ,
  AC_MSG_ERROR("Could not find basic C math library!"))

# Basic compression library
AC_CHECK_LIB([z],[gzwrite],
  ,
  AC_MSG_ERROR("Could not find libz compression library!"))

AC_CHECK_LIB([cpgplot],cpgenv,
             ,
	     AC_MSG_WARN("Could not find libcpgplot plotting library."))

AC_CHECK_LIB([ncurses],[tputs])

AC_CHECK_LIB([readline],[readline],
             ,
             AC_MSG_WARN("Could not find readline line reading library."))

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

AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_LIB([cfitsio],[ffgky],
    ,
    [AC_MSG_ERROR([Couldn't find cfitsio libraries. You will not be able to
                  write 2D confidence images. To compile without cfitsio
                  support, specify --without-cfitsio.])])],
  [AC_MSG_NOTICE([Not compiling CFITSIO image support.])])
AM_CONDITIONAL(CFITSIO,test x$with_cfitsio != xno)

AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_HEADERS([fitsio.h],
  ,
  [AC_MSG_ERROR([Could not validate fitsio.h. Make sure the cfitsio
                 development files are installed and available.])])])


AC_ARG_WITH([debug], 
  [AS_HELP_STRING([--with-debug],
      [Disable optimization/enable debug symbols])],
  [AC_SUBST(OPTFLAGS,"-O0 -g")],
  [with_debug=no])

AM_CONDITIONAL(DEBUG,test x$with_debug = xyes)

AC_CHECK_LIB([gslcblas], [cblas_dnrm2],
  ,
  [AC_MSG_ERROR("Could not find GSL library gslcblas.")])

AC_CHECK_LIB([gsl], [gsl_vector_alloc],
  ,
  [AC_MSG_ERROR("Could not find the main GSL library.")])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt strstr])

AC_CONFIG_SUBDIRS([sciutils])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
