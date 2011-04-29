#
# SYNOPSIS
#
#   AX_GUESS_PATH_HEADER([foo.h])
#
# DESCRIPTION
#
#   Search for header foo.h in -Ipath's found in CPPFLAGS and CFLAGS and set FOO_H_PATH to
#   the full directory path where foo.h was found.
#   If no header is found in the paths given in CPPFLAGS and CFLAGS, then lastly it looks in /usr/local/include.
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

AC_DEFUN([AX_GUESS_PATH_HEADER],
[
    function cw_search_header_path
    {
      n=2
      while test $n -le [$]#; do
	eval arg=\$"$n"
	case "$arg" in
	  -I*)
	    path="`echo "$arg" | sed -e 's/-I//'`"
	    if test -e "$path/$1"; then
	      echo "$path"
	      return
	    fi
	    ;;
	esac
	n=$((n+1))
      done
      if test -e "/usr/local/include/$1"; then
        echo "/usr/local/include"
      fi
    }

    have_realpath=`which realpath`

    cw_headername_uppercase=`echo "m4_toupper([$1])" | sed -e 's/[[^A-Z]]/_/g'`
    AC_CACHE_CHECK([if we can find [$1]], [cw_cv_"$[]cw_headername_uppercase"_path],
    [
      cw_header_path=`eval cw_search_header_path [$1] $CPPFLAGS $CFLAGS`
      if test -n "$cw_header_path"; then
        if test "x$have_realpath" != "x"; then
            eval cw_cv_"$cw_headername_uppercase"_path=`realpath -s "$cw_header_path"`
        else
            eval cw_cv_"$cw_headername_uppercase"_path="$cw_header_path"
        fi
      else
        eval cw_cv_"$cw_headername_uppercase"_path="no"
      fi
    ])
    if eval test \"\$cw_cv_"$cw_headername_uppercase"_path\" = "no"; then
      eval "$cw_headername_uppercase"_PATH=""
    else
      eval "$cw_headername_uppercase"_PATH=\"\$cw_cv_"$cw_headername_uppercase"_path\"
    fi
])
