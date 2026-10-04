# Modified version of AC_OPENMP from autoconf 2.73 c.m4 file (Dima Pasechnik, 2026/03)
#
# Check which options need to be passed to the C compiler to support OpenMP.
# Set the OPENMP_CFLAGS / OPENMP_CXXFLAGS / OPENMP_FFLAGS variable to these
# options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports OpenMP.
#
# For each candidate option, we do a compile test first, then a link test;
# if the compile test succeeds but the link test fails, that means we have
# found the correct option but it doesn't work because the libraries are
# broken.  (This can happen, for instance, with SunPRO C and a bad combination
# of operating system patches.)
#
# Several of the options in our candidate list can be misinterpreted by
# compilers that don't use them to activate OpenMP support; for example,
# many compilers understand "-openmp" to mean "write output to a file
# named 'penmp'" rather than "enable OpenMP".  We can't completely avoid
# the possibility of clobbering files named 'penmp' or 'mp' in configure's
# working directory; therefore, this macro will bomb out if any such file
# already exists when it's invoked.
AC_DEFUN([AX_OPENMP],dnl hackish rename s/AC_/AX_
[AC_REQUIRE([_AC_OPENMP_SAFE_WD])]dnl
[AC_ARG_ENABLE([openmp],
   [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])]dnl
[
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='not found'
      dnl Try these flags:
      dnl   (on by default)      ''
      dnl   GCC >= 4.2           -fopenmp
      dnl   SunPRO C             -xopenmp
      dnl   Intel C              -openmp
      dnl   SGI C, PGI C         -mp
      dnl   Tru64 Compaq C       -omp
      dnl   IBM XL C (AIX, Linux) -qsmp=omp
      dnl   Cray CCE             -homp
      dnl   NEC SX               -Popenmp
      dnl   Lahey Fortran (Linux)  --openmp
      dnl   Apple (Darwin) clang   -Xpreprocessor -fopenmp 
  ax_openmp_flags=":-fopenmp:-Xpreprocessor -fopenmp:-xopenmp:-openmp:-mp:-omp:-qsmp=omp:-homp:-Popenmp:--openmp"
  ac_save_ax_openmp_IFS="$IFS"; IFS=":"
  ax_openmp_old_libs=$LIBS

     for ac_option in $ax_openmp_flags; do
        IFS="$ac_save_ax_openmp_IFS"
        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
        AC_COMPILE_IFELSE([_AC_LANG_OPENMP],
          [
           for ax_openmp_lib in '' -lomp; do
             LIBS="${LIBS} $ax_openmp_lib "
             AC_LINK_IFELSE([_AC_LANG_OPENMP],
               [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option
                ac_cv_[]_AC_LANG_ABBREV[]_omplib=$ax_openmp_lib
                break],
               [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'])
           done
          ])
        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS

        if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'not found'; then
          break
        fi
        IFS=":"
      done

  IFS="$ac_save_ax_openmp_IFS"
  LIBS="${ax_openmp_old_libs}"

      if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" = 'not found'; then
        ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
      elif test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" = ''; then
        ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'
      fi
      dnl _AC_OPENMP_SAFE_WD checked that these files did not exist before we
      dnl started probing for OpenMP support, so if they exist now, they were
      dnl created by the probe loop and it's safe to delete them.
      rm -f penmp mp])
    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'unsupported' && \
       test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != 'none needed'; then
      OPENMP_[]_AC_LANG_PREFIX[]FLAGS="$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp"
      OPENMP_[]_AC_LANG_PREFIX[]LIB="$ac_cv_[]_AC_LANG_ABBREV[]_omplib"
    fi
  fi
  AC_SUBST([OPENMP_]_AC_LANG_PREFIX[FLAGS])
  AC_SUBST([OPENMP_]_AC_LANG_PREFIX[LIB])
])
