#! /bin/sh

# Clueless user check.
if test ! -d .hg -a -f configure; then
  echo "You only need to run './autogen.sh' when you checked out this project using git."
  echo "Just run ./configure [--help]."
  echo "If you insist on running it, then first remove the 'configure' script."
  exit 0
fi

# Installation prefix. This has to match what was used while installing cwautomacros of course.
PREFIX=${CWAUTOMACROSPREFIX-/usr}

if test ! -f $PREFIX/share/cwautomacros/scripts/autogen.sh; then
  echo "$0: $PREFIX/share/cwautomacros/scripts/autogen.sh: No such file or directory"
  echo "$0: This project needs 'cwautomacros'. See http://cwautomacros.berlios.de/"
  exit 126
fi

exec $PREFIX/share/cwautomacros/scripts/autogen.sh
