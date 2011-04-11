#
# SYNOPSIS
#
#   AX_FIND_LIB([foo])
#
# DESCRIPTION
#
#   Search for library foo in -Lpath's found in LDFLAGS and set LIBFOO_PATH to
#   the full directory path where it was found. If no library is found then
#   lastly it looks in /usr/local/lib.
#
# LAST MODIFICATION
#
#   2011-04-10
#
# COPYLEFT
#
#   Copyright (c) 2011 Carlo Wood <carlo@alinoe.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_FIND_LIB],
[
    function search_library_path
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
      if test -e /usr/local/lib/lib$1.so; then
        echo "/usr/local/lib"
      fi
    }

    ax_libname_uppercase="m4_toupper([$1])"
    ax_searchpath=`eval search_library_path [$1] $LDFLAGS`
    if test -n "$ax_searchpath"; then
      eval LIB"$ax_libname_uppercase"_PATH=`realpath -s $ax_searchpath`
    else
      eval LIB"$ax_libname_uppercase"_PATH=""
    fi
])
