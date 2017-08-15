#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(SciUtils, 0.9, jacocluster-discuss@users.sourceforge.net)
AC_CONFIG_SRCDIR([readdata.c])
AM_INIT_AUTOMAKE([gnu])
AC_CONFIG_HEADER([config.h])
AC_PROG_LIBTOOL

OPTFLAGS=-O3
CFLAGS="-Wall -Wfloat-equal -W -Wshadow -Wbad-function-cast -Wredundant-decls -Wdisabled-optimization  -std=c99 -funroll-loops"
LDFLAGS="$LDFLAGS -L/usr/local/lib -L$HOME/lib -L/opt/local/lib"
CPPFLAGS="-I/opt/local/include -I/usr/local/include -I$HOME/include $CPPFLAGS "

AC_ARG_WITH([debug], 
  [AS_HELP_STRING([--with-debug],
      [Disable optimization/enable debug symbols])],
  [AC_SUBST(OPTFLAGS,"-O0 -g")],
  [with_debug=no])

# Checks for programs.
# AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_CC_C_O

# Check required headers
AC_CHECK_HEADERS([gsl/gsl_version.h stdlib.h string.h strings.h unistd.h])


# Basic math library
AC_CHECK_LIB([m],[sqrt],
  ,
  AC_MSG_ERROR("Could not find basic C math library!"))

AC_ARG_WITH([cfitsio], 
  [AS_HELP_STRING([--without-cfitsio],
      [Disable cfitsio (and standalone JACO)])],
  [],
  [with_cfitsio=yes])

LIBCFITSIO= 
AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_LIB([cfitsio],[ffgky],
    [AC_DEFINE([HAVE_CFITSIO],1,[ Define 1 if compiling CFITSIO support ])
     AC_SUBST(LIBCFITSIO,"-lcfitsio")],
    [AC_MSG_WARN([Couldn't find cfitsio libraries. This is not 
                  a problem unless you need the X-ray astronomy routines.])])],
  [])
AM_CONDITIONAL(CFITSIO,test x$with_cfitsio != xno)

AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_HEADERS([fitsio.h],
  ,
  [AC_MSG_ERROR([Could not validate fitsio.h. Make sure the cfitsio
                 development files are installed and available.])])])


AC_CHECK_LIB([gslcblas], [cblas_dnrm2],
  ,
  [AC_MSG_ERROR("Could not find GSL library gslcblas.")])

AC_CHECK_LIB([gsl], [gsl_vector_alloc],
  ,
  [AC_MSG_ERROR("Could not find the main GSL library.")])

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt strstr])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
