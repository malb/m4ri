# ===========================================================================
#             http://autoconf-archive.cryp.to/ax_cache_size.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CACHE_SIZE
#
# DESCRIPTION
#
#   Find L1 and L2 caches size by reading the corresponding file on UNIX or
#   by requesting cpuid. The results are available in the substituted variables
#   M4RI_CPU_L1_CACHE and M4RI_CPU_L2_CACHE.
#
#   This macro depends on AX_GCC_X86_CPUID, AC_PROG_SED, and AX_CPU_VENDOR.
#
# LAST MODIFICATION
#
#   2011-04-11
#
# COPYLEFT
#
#   Copyright (c) 2008 Christophe Tournayre <turn3r@users.sourceforge.net>
#
# Patched by:
#
#   Copyright (c) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#   Copyright (c) 2008 Arnaud Bergeron <abergeron@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CACHE_SIZE],
[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AX_GCC_X86_CPUID])
    AC_REQUIRE([AX_CPU_VENDOR])

    AX_CPU_VENDOR

    ax_l1_size=
    ax_l2_size=

    #Check if the variable is present
    if test -e /sys/devices/system/cpu/cpu0/cache/index0/size; then
        for idx in `seq 0 3`; do
            if test -e /sys/devices/system/cpu/cpu0/cache/index$idx/size ; then
                level=`cat /sys/devices/system/cpu/cpu0/cache/index$idx/level`
                size=`cat /sys/devices/system/cpu/cpu0/cache/index$idx/size`
                eval CPU0\_L$level\_CACHE="$size"
            fi
        done

        ax_l1_size=$CPU0_L1_CACHE
        ax_l2_size=$CPU0_L2_CACHE
        ax_l3_size=$CPU0_L3_CACHE

    else
      if test "x$ax_cv_cpu_vendor" != "xUnknown"; then
        #Or use CPUID
	AX_GCC_X86_CPUID(0x80000000)
	cpu_exthigh=`echo $ax_cv_gcc_x86_cpuid_0x80000000 | cut -d ":" -f 1`
	if test "x$cpu_exthi" > "x80000004"; then
          AX_GCC_X86_CPUID(0x80000005) # For L1 cache
          l1_hexval=`echo $ax_cv_gcc_x86_cpuid_0x80000005 | cut -d ":" -f 4`
          ax_l1_size=$((0x$l1_hexval >> 24))
        fi

	if test "x$cpu_exthi" > "x80000005"; then
          AX_GCC_X86_CPUID(0x80000006) # For L2 cache
          l2_hexval=`echo $ax_cv_gcc_x86_cpuid_0x80000006 | cut -d ":" -f 3`
          ax_l2_size=$((0x$l2_hexval >> 16))
        fi

	if test "x$cpu_exthi" > "x80000005"; then
          AX_GCC_X86_CPUID(0x80000006) # For L3 cache
          l2_hexval=`echo $ax_cv_gcc_x86_cpuid_0x80000006 | cut -d ":" -f 4`
          ax_l2_size=$((0x$l2_hexval >> 18))*512
        fi

      fi

      #Or use sysctl
      sysctl_exe=
      if test -x /usr/sbin/sysctl ; then
	sysctl_exe=/usr/sbin/sysctl
      elif test -x /sbin/sysctl ; then
	sysctl_exe=/sbin/sysctl
      fi
      if test -n "$sysctl_exe"; then
	if test -z "$ax_l2_size" -o "$ax_l2_size" = "0"; then
          sysctl_out=`$sysctl_exe -n hw.l2cachesize 2>/dev/null`;
          if test ! -z "$sysctl_out"; then
            ax_l2_size=$(($sysctl_out / 1024))
          fi;

	fi
	if test -z "$ax_l1_size" -o "$ax_l1_size" = "0" ; then
          sysctl_out=`$sysctl_exe -n hw.l1dcachesize 2>/dev/null`;
          if test ! -z "$sysctl_out"; then
  	    ax_l1_size=$(($sysctl_out / 1024))
          fi;
	fi
	if test -z "$ax_l1_size" -o "ax_l1_size" = "0" ; then
          sysctl_out=`$sysctl_exe -n hw.l1cachesize 2>/dev/null`;
          if test ! -z "$sysctl_out"; then
             ax_l1_size=$(($sysctl_out / 1024))
          fi;
        fi
      fi
    fi

    test -z "$ax_l1_size" && ax_l1_size=0
    test -z "$ax_l2_size" && ax_l2_size=0
    test -z "$ax_l3_size" && ax_l3_size=$ax_l2_size

    # Keep only digits if there is a unit (ie 1024K -> 1024) and convert in Bytes
    AC_MSG_CHECKING(the L1 cache size)
    ax_l1_size=`echo $ax_l1_size | $SED 's/\([[0-9]]\)[[A-Za-z]]$/\1/g'`
    ax_l1_size=$(($ax_l1_size*1024))
    AC_MSG_RESULT( $ax_l1_size Bytes)

    AC_MSG_CHECKING(the L2 cache size)
    ax_l2_size=`echo $ax_l2_size | $SED 's/\([[0-9]]\)[[A-Za-z]]$/\1/g'`
    ax_l2_size=$(($ax_l2_size*1024))
    AC_MSG_RESULT( $ax_l2_size Bytes)

    AC_MSG_CHECKING(the L3 cache size)
    ax_l3_size=`echo $ax_l3_size | $SED 's/\([[0-9]]\)[[A-Za-z]]$/\1/g'`
    ax_l3_size=$(($ax_l3_size*1024))
    AC_MSG_RESULT( $ax_l3_size Bytes)

    M4RI_CPU_L1_CACHE=${ax_l1_size}
    M4RI_CPU_L2_CACHE=${ax_l2_size}
    M4RI_CPU_L3_CACHE=${ax_l3_size}
    AC_SUBST(M4RI_CPU_L1_CACHE)
    AC_SUBST(M4RI_CPU_L2_CACHE)
    AC_SUBST(M4RI_CPU_L3_CACHE)
])
