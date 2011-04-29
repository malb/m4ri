#
# SYNOPSIS
#
#   AX_GUESS_PATH_LIB([foo])
#
# DESCRIPTION
#
#   Search for library foo in -Lpath's found in LDFLAGS and set LIBFOO_PATH to
#   the full directory path where libfoo.so was found.
#   If no library is found in paths given in LDFLAGS, then lastly it looks in /usr/local/lib.
#
# LAST MODIFICATION
#
#   2011-04-11
#
# COPYLEFT
#
#   Copyright (c) 2011 Carlo Wood <carlo@alinoe.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_GUESS_PATH_LIB],
[
    function cw_search_library_path
    {
      n=2
      while test $n -le [$]#; do
	eval arg=\$"$n"
	case "$arg" in
	  -L*)
	    path="`echo "$arg" | sed -e 's/-L//'`"
	    if test -e "$path/lib$1.so"; then
	      echo "$path"
	      return
	    fi
	    ;;
	esac
	n=$((n+1))
      done
      if test -e "/usr/local/lib/lib$1.so"; then
        echo "/usr/local/lib"
      fi
    }

    have_realpath=`which realpath`

    cw_libname_uppercase="m4_toupper([$1])"
    AC_CACHE_CHECK([if we can find lib[$1].so], [cw_cv_lib"$[]cw_libname_uppercase"_path],
    [
      cw_library_path=`eval cw_search_library_path [$1] $LDFLAGS`
      if test -n "$cw_library_path"; then
        if test "x$have_realpath" != "x"; then
           eval cw_cv_lib"$cw_libname_uppercase"_path=`realpath -s "$cw_library_path"`
         else
           eval cw_cv_lib"$cw_libname_uppercase"_path="$cw_library_path"
        fi
      else
        eval cw_cv_lib"$cw_libname_uppercase"_path="no"
      fi
    ])
    if eval test \"\$cw_cv_lib"$cw_libname_uppercase"_path\" = "no"; then
      eval LIB"$cw_libname_uppercase"_PATH=""
    else
      eval LIB"$cw_libname_uppercase"_PATH=\"\$cw_cv_lib"$cw_libname_uppercase"_path\"
    fi
])
