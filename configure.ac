dnl Jasjeet S. Sekhon <sekhon@berkeley.edu> http://sekhon.berkeley.edu
dnl Current version: 2013-10-28
dnl First version: 2006-09-17
dnl Process this file with autoconf to produce a configure script.
dnl autoconf configure.ac > configure

dnl TODO: optimize options for clang, e.g, -O4

dnl Note: to manually switch the compiler and link options one can use
dnl MAKEFLAGS="CC=icc" R CMD SHLIB *.c* amongst other ways.  This 
dnl script handles gcc well.  For the Intel compilers see makefile.icc
dnl in the "extras" directory, and for Sun Studio 12 see makefile.sun.
dnl Modify these makefiles for your environment.

dnl Matching Version 4.8-0+ we now check out flags along with
dnl RC(XX)FLAGS. Ripley noted an issue with -ffast-math and Fedora 16:
dnl -ffast-math and -pedantic (which is an option that R adds) are
dnl incompatible in that OS. Update this script to also check use with RCXXFLAGS
dnl and test script updated with #include <cmath>.

dnl Matching Version 3.6-0+ we switched to configure.ac instead of the
dnl configure hack.

dnl Matching Version 4.8-2: -ffriend-injection is no longer supported in recent versions of gcc
dnl Matching Version 3.4.4+: we need to add the -ffriend-injection
dnl argument for g++ versions 4.1+ because ARM-style name-injection of
dnl friend declarations is no longer the default.  We check for the
dnl flag directly instead of relying on gcc --version information because of 
dnl backporting issues---e.g., Gentoo.

dnl We also setup optimization flags.  This significantly improves performance.

dnl On Windows configure does *not* need to be run.  On Windows,
dnl Makevars.win should simply be used because we can assume that we are
dnl using the mingw32 compiler

dnl -fstrict-aliasing -freorder-blocks -fsched-interblock is added for OS X 
dnl Already included in the Linux gcc -O2 option, but does no harm on Linux

AC_INIT(Matching, 4.0.8+, sekhon@berkeley.edu)

dnl Find the compilers to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
RCFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
RCXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

dnl MYCFLAGS0 should work on most compilers and MYCFLAGS1 should work on most GNU C compilers.  
dnl MYCFLAGS2 may not (generally needs gcc 4.0+).   The last three options in MYCFLAGS2 are the OS X 
dnl addons referenced above.

MYCFLAGS0="-O3"
MYCFLAGS1="-finline-functions -funroll-loops -fexpensive-optimizations"
MYCFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock"

dnl -funsafe-loop-optimizations for gcc 4.2+; the following doesn't work any better in any case
dnl MYCFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock -mfpmath=sse,387 -fgcse -fgcse-after-reload -fgcse-lm -fgcse-sm -fsched-spec-load-dangerous -funsafe-loop-optimizations"



MYCXXFLAGS0="-O3"
MYCXXFLAGS1="-finline-functions -funroll-loops -fexpensive-optimizations"
MYCXXFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock"

dnl -funsafe-loop-optimizations for gcc 4.2+; the following doesn't work any better in any case
dnl MYCXXFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock -mfpmath=sse,387 -fgcse -fgcse-after-reload -fgcse-lm -fgcse-sm -fsched-spec-load-dangerous -funsafe-loop-optimizations"


AC_PROG_CC

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS0])
AC_LANG([C])
CFLAGS="$MYCFLAGS0 $RCFLAGS"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char b[10];]])], [mycflags0_set=yes], [mycflags0_set=no])    

dnl AC_COMPILE_IFELSE([[char b[10];]], [], [AC_MSG_ERROR([you lose])])
dnl AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char b[10];]])], [], [AC_MSG_ERROR([you lose])])

AC_MSG_RESULT($mycflags0_set)    

if test $mycflags0_set = no; then    
   PKG_CLFAGS=""
else
   PKG_CFLAGS=$MYCFLAGS0
fi

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS1])
AC_LANG([C])
CFLAGS="$MYCFLAGS1 $RCFLAGS -pedantic"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char b[10];]])], [mycflags1_set=yes], [mycflags1_set=no])    
AC_MSG_RESULT($mycflags1_set)    

if test $mycflags1_set = yes; then    
   PKG_CFLAGS="$PKG_CFLAGS $MYCFLAGS1"
fi

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS2])
AC_LANG([C])
CFLAGS="$MYCFLAGS2 $RCFLAGS"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[char b[10];]])], [mycflags2_set=yes], [mycflags2_set=no])    
AC_MSG_RESULT($mycflags2_set)    

if test $mycflags2_set = yes; then    
   PKG_CFLAGS="$PKG_CFLAGS $MYCFLAGS2"
fi

AC_PROG_CXX

dnl if test $GXX = no; then
dnl    AC_MSG_WARN([**** GNU g++ IS HIGHLY RECOMMENDED.  I will continue, but this may not work. ****])
dnl fi

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS0])
AC_LANG([C++])
CXXFLAGS="$MYCXXFLAGS0 $RCXXFLAGS"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <cmath>
                    int main() {
                      std::sqrt(5.1);
                      return 0;
                    }]])], [mycxxflags0_set=yes], [mycxxflags0_set=no])    
AC_MSG_RESULT($mycxxflags0_set)    

if test $mycxxflags0_set = no; then    
  PKG_CXXLFAGS=""
else
  PKG_CXXFLAGS=$MYCXXFLAGS0
fi

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS1])
AC_LANG([C++])
CXXFLAGS="$MYCXXFLAGS1 $RCXXFLAGS -pedantic"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <cmath>
                    int main() {
                      std::sqrt(5.1);
                      return 0;
                    }]])], [mycxxflags1_set=yes], [mycxxflags1_set=no])    
AC_MSG_RESULT($mycxxflags1_set)    

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS2])
AC_LANG([C++])
CXXFLAGS="$MYCXXFLAGS2 $RCXXFLAGS"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <cmath>
                    int main() {
                      std::sqrt(5.1);
                      return 0;
                    }]])], [mycxxflags2_set=yes], [mycxxflags2_set=no])    
AC_MSG_RESULT($mycxxflags2_set)    

if test $mycxxflags1_set = yes; then    
   PKG_CXXFLAGS="$PKG_CXXFLAGS $MYCXXFLAGS1"
fi

if test $mycxxflags2_set = yes; then    
   PKG_CXXFLAGS="$PKG_CXXFLAGS $MYCXXFLAGS2"
fi

dnl -ffrend-injection is no longer supported
dnl AC_MSG_CHECKING([whether $CXX accepts -ffriend-injection])
dnl AC_LANG([C++])
dnl CXXFLAGS="-ffriend-injection $RCXXFLAGS"
dnl AC_COMPILE_IFELSE([[char b[10];]], [ffriend_set=yes], [ffriend_set=no])
dnl AC_MSG_RESULT($ffriend_set)

dnl if test $ffriend_set = yes; then    
dnl    PKG_CXXFLAGS="$PKG_CXXFLAGS -ffriend-injection"
dnl fi

dnl Are we on OSX?
dnl OS X specific modifications are no longer needed: Oct 28, 2013
dnl If so: use dmalloc (NO LONGER DONE) and don't compile our own cblas headers (still done).
dnl AC_MSG_CHECKING([whether we are on Darwin])
dnl UNAME=`uname -a`
dnl echo $UNAME | grep -i Darwin > /dev/null 2>&1
dnl if test "$?" = 0 ; then
dnl  AC_MSG_RESULT(yes)
dnl cp inst/extras/malloc.c src/malloc.c
dnl  mv src/cblas_dgemm.c inst/extras
dnl  mv src/cblas.h inst/extras
dnl  mv src/cblas_dasum.c inst/extras
dnl  mv src/cblas_dscal.c inst/extras
dnl else
dnl  AC_MSG_RESULT(no)
dnl fi

AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_CXXFLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

