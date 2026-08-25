#!/bin/sh

(
  echo 1 /proc/sys/kernel/perf_user_access
  echo 2 /proc/sys/kernel/perf_event_paranoid
  echo 2 /sys/devices/cpu/rdpmc
  echo 2 /sys/devices/cpu_core/rdpmc
  echo 2 /sys/devices/cpu_atom/rdpmc
  echo 0 /proc/sys/kernel/nmi_watchdog
) | while read x fn
do
  [ -e $fn ] && echo $x > $fn
done
