dnl Process this file with autoconf to produce a configure script, or use autoreconf -i

CFLAGS="$CFLAGS"

AC_INIT(itcsimlib, version-0.0)
AC_PROG_CC
AC_LANG(C)

AC_CHECK_HEADERS(stdlib)
AC_CHECK_HEADERS(math)

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_SUBST(LIBS)
AC_SUBST(CFLAGS)

AC_OUTPUT(makefile)