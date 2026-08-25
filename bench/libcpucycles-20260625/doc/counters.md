Currently libcpucycles supports the following cycle counters. Some
cycle counters are actually other forms of counters that libcpucycles
scales to imitate a cycle counter. There is
[separate documentation](selection.html)
for how libcpucycles makes a choice of cycle counter.
There is a `cpucycles-open` command that, when run as `root`, often
succeeds in enabling higher-quality cycle counters until reboot. See the
[security considerations](security.html) regarding enabling or disabling
counters and regarding overclocking (Turbo Boost, Turbo Core, etc.).

`amd64-perfpmcff`: Requires a 64-bit Intel platform; the Linux
`perf_event` interface; `/proc/sys/kernel/perf_event_paranoid` at most 2;
any existing files among `/sys/devices/{cpu,cpu_core,cpu_atom}/rdpmc` at
least 2. Accesses a cycle counter through RDPMC (using Intel-specific
counter `0x40000001`, which is why this ends up not running on AMD CPUs).
This counter runs at the clock frequency of the CPU core.
Starting with libcpucycles version 20260625, this counter compensates
for a `perf_event` bug in handling Intel E-cores.

`amd64-perfpmc`: Requires a 64-bit AMD platform; the Linux `perf_event`
interface; `/proc/sys/kernel/perf_event_paranoid` at most 2; any
existing files among `/sys/devices/{cpu,cpu_core,cpu_atom}/rdpmc` at
least 2. Accesses a cycle counter through RDPMC (using the `perf_event`
index); to avoid various `amd64-perfpmcff` complications, actively
refuses to run on Intel. This counter runs at the clock frequency of the
CPU core.

`amd64-pmcff`: Requires a 64-bit Intel platform and an OS that enables
user-level RDPMC access. Accesses a cycle counter through RDPMC (using
Intel-specific counter `0x40000001`, which is why this ends up not running
on AMD CPUs). This counter runs at the clock frequency of the CPU core.

`amd64-tsc`, `amd64-tscasm`: Requires a 64-bit Intel/AMD platform.
Requires RDTSC to be enabled, which it is by default. Uses RDTSC to
access the CPU's time-stamp counter. On current CPUs, this is an
off-core clock rather than a cycle counter, but it is typically a very
fast off-core clock, making it adequate for seeing cycle counts if
overclocking and underclocking are disabled. The difference between
`tsc` and `tscasm` is that `tsc` uses the compiler's `__rdtsc()` while
`tscasm` uses inline assembly.

`arm32-cortex`: Requires a 32-bit ARMv7-A platform. Uses
`mrc p15, 0, %0, c9, c13, 0` to read the cycle counter. Requires user
access to the cycle counter, which is not enabled by default but can (sometimes) be
enabled under Linux via
[a kernel module](https://github.com/thoughtpolice/enable_arm_pmu).
This counter is natively 32 bits, but libcpucycles watches how the
counter and `gettimeofday` increase to compute a 64-bit extension of the
counter.

`arm32-1176`: Requires a 32-bit ARM1176 platform. Uses
`mrc p15, 0, %0, c15, c12, 1` to read the cycle counter. Requires user
access to the cycle counter, which is not enabled by default but can (sometimes) be
enabled under Linux via
[a kernel module](https://bench.cr.yp.to/cpucycles/n810.html).
This counter is natively 32 bits, but libcpucycles watches how the
counter and `gettimeofday` increase to compute a 64-bit extension of the
counter.

`arm64-perfpmc`: Requires a 64-bit ARMv8-A platform; the Linux
`perf_event` interface; on current Linux kernels, 1 in
`/proc/sys/kernel/perf_user_access`; on older Linux kernels,
[a kernel module](https://github.com/rdolbeau/enable_arm_pmu)
to enable PMCCNTR access. Uses `mrs %0, PMCCNTR_EL0` to read the cycle
counter.

`arm64-pmc`: Requires a 64-bit ARMv8-A platform and an OS that enables
user-level PMCCNTR access, for example an older Linux kernel with
[a kernel module](https://github.com/rdolbeau/enable_arm_pmu).
Uses `mrs %0, PMCCNTR_EL0` to read the cycle counter.

`arm64-vct`: Requires a 64-bit ARMv8-A platform. Uses
`mrs %0, CNTVCT_EL0` to read a "virtual count" timer. This is an
off-core clock, typically running at 24MHz. Results are scaled by
libcpucycles.

`loong64-rdtime`: Requires a 64-bit LoongArch platform.
Uses `rdtime.d` to read the "constant frequency timer".
This is an off-core clock, typically running at 100MHz.

`mips64-cc`: Requires a 64-bit MIPS platform. (Maybe the same code would
also work as `mips32-cc`, but this has not been tested yet.) Uses RDHWR
to read the hardware cycle counter (hardware register 2 times a constant
scale factor in hardware register 3). This counter is natively 32 bits,
but libcpucycles watches how the counter and `gettimeofday` increase to
compute a 64-bit extension of the counter.

`ppc32-mftb`: Requires a 32-bit PowerPC platform. Uses `mftb` and
`mftbu` to read the "time base". This is an off-core clock, typically
running at 24MHz.

`ppc64-mftb`: Requires a 64-bit PowerPC platform. Uses `mftb` and
`mftbu` to read the "time base". This is an off-core clock, typically
running at 24MHz.

`riscv32-perfrdcycle`: Requires a 32-bit RISC-V platform; the Linux
`perf_event` interface; on current Linux kernels, 1 in
`/proc/sys/kernel/perf_user_access`. Uses `rdcycle` and `rdcycleh` to
read a cycle counter.

`riscv32-rdcycle`: Requires a 32-bit RISC-V platform and an OS that
enables user-level `rdcycle` access. Uses `rdcycle` and `rdcycleh` to
read a cycle counter.

`riscv64-perfrdcycle`: Requires a 64-bit RISC-V platform; the Linux
`perf_event` interface; on current Linux kernels, 1 in
`/proc/sys/kernel/perf_user_access`. Uses `rdcycle` to read a cycle
counter.

`riscv64-rdcycle`: Requires a 64-bit RISC-V platform and an OS that
enables user-level `rdcycle` access. Uses `rdcycle` to read a cycle
counter.

`s390x-stckf`: Requires a 64-bit z/Architecture platform. Uses `stckf`
to read the TOD clock, which is documented to run at 4096MHz. On the
z15, this looks like a doubling of an off-core 2048MHz clock. Results
are scaled by libcpucycles.

`sparc64-rdtick`: Requires a 64-bit SPARC platform. Uses `rd %tick`
to read a cycle counter.

`x86-tsc`, `x86-tscasm`: Same as `amd64-tsc` and `amd64-tscasm`, but
for 32-bit Intel/AMD platforms instead of 64-bit Intel/AMD platforms.

`default-gettimeofday`: Reasonably portable. Resolution is limited to 1
microsecond. Results are scaled by libcpucycles.

`default-mach`: Requires an OS with `mach_absolute_time()`. Typically
runs at 24MHz. Results are scaled by libcpucycles.

`default-monotonic`: Requires `CLOCK_MONOTONIC`. Reasonably portable,
although might fail on older systems where `default-gettimeofday` works.
Resolution is limited to 1 nanosecond. Can be almost as good as a cycle
counter, or orders of magnitude worse, depending on the OS and CPU.
Results are scaled by libcpucycles.

`default-perfevent`: Requires the Linux `perf_event` interface, and a
CPU where `perf_event` supports `PERF_COUNT_HW_CPU_CYCLES`. Similar
variations in quality to `default-monotonic`, without the 1-nanosecond
limitation. Starting with libcpucycles version 20260625, this counter
compensates for a `perf_event` bug in handling Intel E-cores.

`default-zero`: The horrifying last resort if nothing else works.

## Examples

These are examples of `cpucycles-info` output on various machines. The
machines named `cfarm*` are from the
[cfarm compile farm](https://portal.cfarm.net/).

A `median` line saying, e.g., `47 +47+28+0+2-5+0+2-5...` means that the
differences between adjacent cycle counts were 47+47, 47+28, 47+0, 47+2,
47−5, 47+0, 47+2, 47−5, etc., with median difference 47. The first few
differences are typically larger because of cache effects.
Current versions of libcpucycles use `iqm` (interquartile mean)
instead of `median`.

`berry0`,
Broadcom BCM2835:
```
cpucycles version 20240114
cpucycles tracesetup 0 arm32-cortex precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 arm32-1176 precision 22 scaling 1.000000 only32 1
cpucycles tracesetup 2 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 1199 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 1200 scaling 1000.000000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1000000000
cpucycles implementation arm32-1176
cpucycles median 720 +942+124+1+0+2+0+0+0+0+0+0+0+0+0+0+0+0+0+0+1+2+0+0+2+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+222+300+1+0+0+2+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 798307692...2045181819 with 1024 loops 12 microseconds
cpucycles observed persecond 915478260...1260523810 with 2048 loops 22 microseconds
cpucycles observed persecond 947809523...1106100000 with 4096 loops 41 microseconds
cpucycles observed persecond 966353658...1129037500 with 8192 loops 81 microseconds
cpucycles observed persecond 988490566...1030019109 with 16384 loops 158 microseconds
cpucycles observed persecond 995169327...1002034063 with 32768 loops 2379 microseconds
cpucycles observed persecond 996871019...1012568691 with 65536 loops 627 microseconds
cpucycles observed persecond 997832134...1004212170 with 131072 loops 1250 microseconds
cpucycles observed persecond 997740918...1000887780 with 262144 loops 5009 microseconds
cpucycles observed persecond 998528349...1001961164 with 524288 loops 5537 microseconds
cpucycles observed persecond 999202882...1001166794 with 1048576 loops 10547 microseconds
```

`berry2`,
Broadcom BCM2836:
```
cpucycles version 20251226
cpucycles tracesetup 0 arm32-cortex precision 310 scaling 1.000000 only32 1
cpucycles tracesetup 1 arm32-1176 precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 254 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 814 scaling 0.900000 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 1200 scaling 900.000000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 900000000
cpucycles implementation default-perfevent
cpucycles iqm 154 +208+9+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+173+172+29+0-4+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 218360000...420565218 with 1024 loops 24 microseconds
cpucycles observed persecond 327781250...402566667 with 2048 loops 31 microseconds
cpucycles observed persecond 423204081...468085107 with 4096 loops 48 microseconds
cpucycles observed persecond 490369047...517231708 with 8192 loops 83 microseconds
cpucycles observed persecond 544033112...560208054 with 16384 loops 150 microseconds
cpucycles observed persecond 450227397...457187328 with 32768 loops 364 microseconds
cpucycles observed persecond 584518716...588967800 with 65536 loops 560 microseconds
cpucycles observed persecond 592214995...594508598 with 131072 loops 1106 microseconds
cpucycles observed persecond 595894545...597049591 with 262144 loops 2199 microseconds
cpucycles observed persecond 598008895...598558878 with 524288 loops 4383 microseconds
cpucycles observed persecond 595624105...595918873 with 1048576 loops 8802 microseconds
```

`pi3aplus`,
Broadcom BCM2837B0:
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 9 scaling 1.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1083 scaling 73.000000 only32 0
cpucycles tracesetup 3 default-perfevent precision 130 scaling 1.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1102 scaling 1.400000 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 2440 scaling 1400.000000 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1400000000
cpucycles implementation arm64-perfpmc
cpucycles iqm 10 +15+8+8+3+0+0+0+0+0+0+0+7-2-2-2-1+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 1030500000...2158500000 with 2048 loops 3 microseconds
cpucycles observed persecond 1177285714...1688600000 with 4096 loops 6 microseconds
cpucycles observed persecond 1369333333...1662000000 with 8192 loops 11 microseconds
cpucycles observed persecond 1367000000...1499454546 with 16384 loops 23 microseconds
cpucycles observed persecond 1366166666...1429565218 with 32768 loops 47 microseconds
cpucycles observed persecond 1394808510...1427086957 with 65536 loops 93 microseconds
cpucycles observed persecond 1165577777...1177856503 with 131072 loops 224 microseconds
cpucycles observed persecond 1379878947...1387708995 with 262144 loops 379 microseconds
cpucycles observed persecond 1390775862...1394730054 with 524288 loops 753 microseconds
cpucycles observed persecond 1364528952...1366583714 with 1048576 loops 1536 microseconds
```

`pi4b`,
Broadcom BCM2711,
Ubuntu 24.04,
Linux kernel 6.8.0-1036-raspi,
running in power-saving mode (`echo powersave > /sys/devices/system/cpu/cpufreq/policy0/scaling_governor`):
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 9 scaling 1.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1037 scaling 27.750000 only32 0
cpucycles tracesetup 3 default-perfevent precision 162 scaling 1.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1252 scaling 1.500000 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 2540 scaling 1500.000000 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1500000000
cpucycles implementation arm64-perfpmc
cpucycles iqm 9 +0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 423800000...822333334 with 1024 loops 4 microseconds
cpucycles observed persecond 522000000...729166667 with 2048 loops 7 microseconds
cpucycles observed persecond 551466666...654461539 with 4096 loops 14 microseconds
cpucycles observed persecond 586642857...641384616 with 8192 loops 27 microseconds
cpucycles observed persecond 585892857...611981482 with 16384 loops 55 microseconds
cpucycles observed persecond 596163636...609453704 with 32768 loops 109 microseconds
cpucycles observed persecond 598694063...605235024 with 65536 loops 218 microseconds
cpucycles observed persecond 564165591...567133910 with 131072 loops 464 microseconds
cpucycles observed persecond 587983183...589557304 with 262144 loops 891 microseconds
cpucycles observed persecond 585639865...586455562 with 524288 loops 1790 microseconds
cpucycles observed persecond 593200509...593615450 with 1048576 loops 3535 microseconds
```

`pi5`,
Broadcom BCM2712,
Debian 12,
Linux kernel 6.12.25+rpt-rpi-2712,
running in power-saving mode (`echo powersave > /sys/devices/system/cpu/cpufreq/policy0/scaling_governor`):
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 2 scaling 1.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1054 scaling 44.500000 only32 0
cpucycles tracesetup 3 default-perfevent precision 97 scaling 1.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1162 scaling 2.400000 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 3440 scaling 2400.000000 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2400000000
cpucycles implementation arm64-perfpmc
cpucycles iqm 7 +11+12+1+12+1+12+1+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+3-3+2+1-3+2-2+0+0+0+0+0+3-3+2-2+0+0+0+0+3-3+2-2+0+0+0+0+3-3+2-2+0+0+0+0+0+0+0
cpucycles observed persecond 1034000000...2137000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1178428571...1673800000 with 8192 loops 6 microseconds
cpucycles observed persecond 1370416666...1656500000 with 16384 loops 11 microseconds
cpucycles observed persecond 1426434782...1568000000 with 32768 loops 22 microseconds
cpucycles observed persecond 1457222222...1527790698 with 65536 loops 44 microseconds
cpucycles observed persecond 1489897727...1525930233 with 131072 loops 87 microseconds
cpucycles observed persecond 1489687500...1507494253 with 262144 loops 175 microseconds
cpucycles observed persecond 1493823361...1502727794 with 524288 loops 350 microseconds
cpucycles observed persecond 1495885877...1500337626 with 1048576 loops 700 microseconds
```

`bblack`,
TI Sitara XAM3359AZCZ100:
```
cpucycles version 20260625
cpucycles tracesetup 0 arm32-cortex precision 308 scaling 1.000000 only32 1
cpucycles tracesetup 1 arm32-1176 precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 2134 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 2020 scaling 1000.000000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1000000000
cpucycles implementation arm32-cortex
cpucycles iqm 1260 +166+5+19+0+10+0+0+0+0+0+0+0+0+0+0+0+0+10+0+0+0+2+8+0+0+0+0+0+0+2+8+2+8+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+12+8+0+0+10+0+0+0+0+10+0+2+8+0+0
cpucycles observed persecond 625636363...5378111112 with 1024 loops 10 microseconds
cpucycles observed persecond 768600000...1400692308 with 2048 loops 14 microseconds
cpucycles observed persecond 869240000...1162217392 with 4096 loops 24 microseconds
cpucycles observed persecond 938333333...1097976745 with 8192 loops 44 microseconds
cpucycles observed persecond 956057471...1037282353 with 16384 loops 86 microseconds
cpucycles observed persecond 977035502...1018814372 with 32768 loops 168 microseconds
cpucycles observed persecond 987789789...1008755288 with 65536 loops 332 microseconds
cpucycles observed persecond 994892424...1005542554 with 131072 loops 659 microseconds
cpucycles observed persecond 997695057...1002990861 with 262144 loops 1314 microseconds
cpucycles observed persecond 998862529...1000671667 with 524288 loops 6735 microseconds
cpucycles observed persecond 999072352...1000502726 with 1048576 loops 9356 microseconds
```

`titan0`,
Intel Xeon E3-1275 V3:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 37 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1034 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1044 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 112 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1106 scaling 3.500000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 4550 scaling 3500.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3500000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 41 +43+20+20+20-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0-1+2-2+0
cpucycles observed persecond 2064000000...4222500000 with 8192 loops 3 microseconds
cpucycles observed persecond 2752500000...4175250000 with 16384 loops 5 microseconds
cpucycles observed persecond 3283900000...4127875000 with 32768 loops 9 microseconds
cpucycles observed persecond 3158724137...3402185186 with 65536 loops 28 microseconds
cpucycles observed persecond 3294725000...3473026316 with 131072 loops 39 microseconds
cpucycles observed persecond 3450223684...3546000000 with 262144 loops 75 microseconds
cpucycles observed persecond 3449710526...3496933334 with 524288 loops 151 microseconds
cpucycles observed persecond 3483853820...3507775920 with 1048576 loops 300 microseconds
```

`nucnuc`,
Intel Pentium N3700:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 18 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1030 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1040 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 371 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1160 scaling 1.600000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 2650 scaling 1600.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1600000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 58 +15+0+0+0+0+15+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 1064250000...2331500000 with 2048 loops 3 microseconds
cpucycles observed persecond 1389000000...2181000000 with 4096 loops 5 microseconds
cpucycles observed persecond 1501181818...1875777778 with 8192 loops 10 microseconds
cpucycles observed persecond 1495318181...1663350000 with 16384 loops 21 microseconds
cpucycles observed persecond 1563476190...1650600000 with 32768 loops 41 microseconds
cpucycles observed persecond 1561916666...1604512196 with 65536 loops 83 microseconds
cpucycles observed persecond 1589533333...1611306749 with 131072 loops 164 microseconds
cpucycles observed persecond 1492115819...1502056819 with 262144 loops 353 microseconds
cpucycles observed persecond 1596178082...1601633588 with 524288 loops 656 microseconds
cpucycles observed persecond 1598525152...1601257252 with 1048576 loops 1311 microseconds
```

`saber214`,
AMD FX-8350,
running under Xen:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1070 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1080 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1188 scaling 4.013176 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 5043 scaling 4013.176000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 4013176000
cpucycles implementation amd64-tsc
cpucycles iqm 71 +46-1+0+0+1+0+0+0+0+1+0+0+0+0+1+0+0+0+1+0+0+0+0+1+0+0+0+1+0+0+0+0+1+0+0+0+1+0+0+0+0+1+0+0+0+1+0+0+0+0+1+0+0+0+0+1+0+0+0+1+0+0+0+0
cpucycles observed persecond 3416000000...7044000000 with 4096 loops 3 microseconds
cpucycles observed persecond 3442142857...4903000000 with 8192 loops 6 microseconds
cpucycles observed persecond 3697538461...4408181819 with 16384 loops 12 microseconds
cpucycles observed persecond 3840800000...4192913044 with 32768 loops 24 microseconds
cpucycles observed persecond 3916857142...4092553192 with 65536 loops 48 microseconds
cpucycles observed persecond 3956072164...4043757895 with 131072 loops 96 microseconds
cpucycles observed persecond 3996682291...4040984211 with 262144 loops 191 microseconds
cpucycles observed persecond 4006827676...4028971129 with 524288 loops 382 microseconds
cpucycles observed persecond 4006678851...4017726440 with 1048576 loops 765 microseconds
```

`phoenix`,
AMD Ryzen 5 7640HS:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 33 scaling 1.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1011 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1021 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 176 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1126 scaling 4.301000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 5351 scaling 4301.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 4301000000
cpucycles implementation amd64-perfpmc
cpucycles iqm 34 +35+0+0+0+0+0+0+481+0+0+0+0+0+0+0+0+0+0+0+559-1+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 3289200000...5535333334 with 16384 loops 4 microseconds
cpucycles observed persecond 3647777777...4715285715 with 32768 loops 8 microseconds
cpucycles observed persecond 3858705882...4383866667 with 65536 loops 16 microseconds
cpucycles observed persecond 4097937500...4377033334 with 131072 loops 31 microseconds
cpucycles observed persecond 4096968750...4231983871 with 262144 loops 63 microseconds
cpucycles observed persecond 4096484375...4162777778 with 524288 loops 127 microseconds
cpucycles observed persecond 4228379032...4263406505 with 1048576 loops 247 microseconds
```

`freshwrap`,
Intel Core 5 210H,
P-core:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 14 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1049 scaling 0.818452 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1059 scaling 0.818452 only32 0
cpucycles tracesetup 5 default-perfevent precision 142 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1103 scaling 2.200000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 3250 scaling 2200.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2200000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 14 +55+27-3-1+0+3-1+1+0+3+5-3+5-1+5-2+3-1+2-1+0+0-1+0+0+1-1+1+0-1+0+0+1-1+1+0-1+0+0+1-1+1+0-1+0+0+1-1+1+0-1+0+0+1-1+1+0-1+0+0+1-1+1+0
cpucycles observed persecond 1828333333...2773500000 with 16384 loops 5 microseconds
cpucycles observed persecond 1991818181...2450777778 with 32768 loops 10 microseconds
cpucycles observed persecond 2082761904...2308368422 with 65536 loops 20 microseconds
cpucycles observed persecond 2132878048...2245358975 with 131072 loops 40 microseconds
cpucycles observed persecond 2165755813...2218761905 with 262144 loops 85 microseconds
cpucycles observed persecond 2185543750...2213962026 with 524288 loops 159 microseconds
cpucycles observed persecond 2189342679...2203442007 with 1048576 loops 320 microseconds
```

`freshwrap`,
Intel Core 5 210H,
E-core:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 14 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1049 scaling 0.818452 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1059 scaling 0.818452 only32 0
cpucycles tracesetup 5 default-perfevent precision 121 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1143 scaling 2.200000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 3250 scaling 2200.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2200000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 17 +21+22+1+3-3+1-1+2-1-1+1-1+3-3+1-1+3-1-1-1+1-1+1-1+1-1+1-1+2-1-1+2-1-1+1-1+1-1+1-1+1-1+1-1+3-3+2-1-1+1-1+3-3+2-1-1+1-1+1-1+2-1-1+1
cpucycles observed persecond 1060750000...2189000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1275571428...1819800000 with 8192 loops 6 microseconds
cpucycles observed persecond 1475250000...1787300000 with 16384 loops 11 microseconds
cpucycles observed persecond 1525217391...1678571429 with 32768 loops 22 microseconds
cpucycles observed persecond 1541955555...1617651163 with 65536 loops 44 microseconds
cpucycles observed persecond 1575720930...1615285715 with 131072 loops 85 microseconds
cpucycles observed persecond 1581636904...1601716868 with 262144 loops 167 microseconds
cpucycles observed persecond 1593329305...1603534955 with 524288 loops 330 microseconds
cpucycles observed persecond 1595613050...1600729072 with 1048576 loops 658 microseconds
```

`meteor`,
Intel Core Ultra 5 125H,
P-core:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 8 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1033 scaling 0.400641 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1043 scaling 0.400641 only32 0
cpucycles tracesetup 5 default-perfevent precision 143 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1104 scaling 1.200000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 2250 scaling 1200.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1200000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 20 +55+25+0+3-4+3+2+9-7-3-3-3-4+0-2+0-4-1-3-3-4-1-1-6-3-3+2+1+6-2-1+10-12+6+2+4-9+2+0-1-1+40-7+3+0+2+0+4-4-1-6-3-3-6+4+3+5-4+12+1+8+10+12+5
cpucycles observed persecond 938750000...1949000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1055500000...1458500000 with 8192 loops 7 microseconds
cpucycles observed persecond 1057875000...1219928572 with 16384 loops 15 microseconds
cpucycles observed persecond 1176060606...1258387097 with 32768 loops 32 microseconds
cpucycles observed persecond 1169929824...1215672728 with 65536 loops 56 microseconds
cpucycles observed persecond 1190088888...1218943182 with 131072 loops 89 microseconds
cpucycles observed persecond 1189457317...1205123457 with 262144 loops 163 microseconds
cpucycles observed persecond 1193379870...1201712419 with 524288 loops 307 microseconds
cpucycles observed persecond 1196813021...1201098828 with 1048576 loops 598 microseconds
```

`meteor`,
Intel Core Ultra 5 125H,
E-core:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 14 scaling 1.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1059 scaling 0.400641 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1068 scaling 0.400641 only32 0
cpucycles tracesetup 5 default-perfevent precision 124 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1167 scaling 1.200000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 2250 scaling 1200.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1200000000
cpucycles implementation amd64-perfpmcff
cpucycles iqm 17 +1-1+1-1+1-1+1-1+1-1+2-1-1+1-1+1-1+2-1-1+1-1+3-3+1-1+1-1+1-1+2-1-1+1-1+1-1+1-1+1-1+2-1-1+1-1+2-1-1+1-1+2-1-1+1-1+1-1+1-1+1-1+1-1
cpucycles observed persecond 439200000...770666667 with 2048 loops 4 microseconds
cpucycles observed persecond 608857142...878200000 with 4096 loops 6 microseconds
cpucycles observed persecond 638000000...764363637 with 8192 loops 12 microseconds
cpucycles observed persecond 687750000...756136364 with 16384 loops 23 microseconds
cpucycles observed persecond 670959183...701978724 with 32768 loops 48 microseconds
cpucycles observed persecond 691105263...707182796 with 65536 loops 94 microseconds
cpucycles observed persecond 694132275...702160428 with 131072 loops 188 microseconds
cpucycles observed persecond 667098984...671061225 with 262144 loops 393 microseconds
cpucycles observed persecond 697416223...699448000 with 524288 loops 751 microseconds
cpucycles observed persecond 697739188...698754164 with 1048576 loops 1502 microseconds
```

`cfarm14`,
Intel Xeon E5-2620 v3,
Debian 12,
Linux kernel 6.1.0-26-amd64:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1018 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1028 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 84 scaling 1.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1110 scaling 3.200000 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 4240 scaling 3200.000000 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3200000000
cpucycles implementation default-perfevent
cpucycles iqm 74 +50+0+0+0+0+0+0+0+0-2+1+0+2+0+0+0+1+0+0+0+0+0+0-2+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+1+0+1+1+0-2+0+0+2+0-2-2-2
cpucycles observed persecond 541000000...1232500000 with 2048 loops 3 microseconds
cpucycles observed persecond 1053750000...2229000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1385500000...2138750000 with 8192 loops 5 microseconds
cpucycles observed persecond 1833555555...2392857143 with 16384 loops 8 microseconds
cpucycles observed persecond 2192533333...2549230770 with 32768 loops 14 microseconds
cpucycles observed persecond 2344750000...2534730770 with 65536 loops 27 microseconds
cpucycles observed persecond 2475471698...2577333334 with 131072 loops 52 microseconds
cpucycles observed persecond 2522778846...2574598040 with 262144 loops 103 microseconds
cpucycles observed persecond 2507608832...2525190477 with 524288 loops 316 microseconds
cpucycles observed persecond 2583098522...2596502476 with 1048576 loops 405 microseconds
```

`cfarm26`,
Intel Core i5-4570 in 32-bit mode under KVM,
Debian 12.12,
Linux kernel 6.1.0-41-686-pae:
```
cpucycles version 20260625
cpucycles tracesetup 0 x86-tsc precision 1018 scaling 1.000000 only32 0
cpucycles tracesetup 1 x86-tscasm precision 1028 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 731 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 3220 scaling 3.192606 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 4232 scaling 3192.606000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3192606000
cpucycles implementation default-perfevent
cpucycles iqm 779 +35-8+9+15+10-3-7-8-15+10+0-8+4+3-1-1+11+11-1+6+4-3-15-7+3-3+0-7+5+3-8+10+4+3+0-15-3-3-1-14-3+3+10-1-3-3-1-7+3-4+0+0-10+3-8+0+4-1+0-7+4+3-8+0
cpucycles observed persecond 283428571...1430400000 with 1024 loops 6 microseconds
cpucycles observed persecond 417142857...1240800000 with 2048 loops 6 microseconds
cpucycles observed persecond 715571428...1662000000 with 4096 loops 6 microseconds
cpucycles observed persecond 1132500000...2059666667 with 8192 loops 7 microseconds
cpucycles observed persecond 1568000000...2281333334 with 16384 loops 10 microseconds
cpucycles observed persecond 2254933333...2854769231 with 32768 loops 14 microseconds
cpucycles observed persecond 2658360000...3034000000 with 65536 loops 24 microseconds
cpucycles observed persecond 2174032786...2310118645 with 131072 loops 60 microseconds
cpucycles observed persecond 3329126582...3459415585 with 262144 loops 78 microseconds
cpucycles observed persecond 3455059210...3523140000 with 524288 loops 151 microseconds
cpucycles observed persecond 3454417763...3489394040 with 1048576 loops 303 microseconds
```

`cfarm27`,
Intel Core i5-4570 in 32-bit mode under KVM,
Alpine 3.19.9,
Linux kernel 6.6.117-0-lts:
```
cpucycles version 20260625
cpucycles tracesetup 0 x86-tsc precision 1018 scaling 1.000000 only32 0
cpucycles tracesetup 1 x86-tscasm precision 1028 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 608 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 1601 scaling 3.192606 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 4232 scaling 3192.606000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3192606000
cpucycles implementation default-perfevent
cpucycles iqm 606 +72+2-9-3+2+6-4-2+2+2-6-3+5+2-8-2+2+2-6-3+2+2-8-2+2+2-6-3+2+2-8-2+2+2-6-3+2+2-8-2+2+2-6-3+2+2-8-2+2+2-6-3+2+2-8-2+2+2-6-3+2+2-8-2
cpucycles observed persecond 356800000...1256666667 with 1024 loops 4 microseconds
cpucycles observed persecond 537400000...1404666667 with 2048 loops 4 microseconds
cpucycles observed persecond 787166666...1562500000 with 4096 loops 5 microseconds
cpucycles observed persecond 1261142857...2072000000 with 8192 loops 6 microseconds
cpucycles observed persecond 1890444444...2648714286 with 16384 loops 8 microseconds
cpucycles observed persecond 2392285714...2934250000 with 32768 loops 13 microseconds
cpucycles observed persecond 2879304347...3228952381 with 65536 loops 22 microseconds
cpucycles observed persecond 2723040816...2872978724 with 131072 loops 48 microseconds
cpucycles observed persecond 3326113924...3432428572 with 262144 loops 78 microseconds
cpucycles observed persecond 3386496774...3440797386 with 524288 loops 154 microseconds
cpucycles observed persecond 3202069696...3226439025 with 1048576 loops 329 microseconds
```

`cfarm29`,
IBM POWER9,
Debian 13.1,
Linux kernel 6.12.57+deb13-powerpc64le:
```
cpucycles version 20260105
cpucycles tracesetup 0 ppc64-mftb precision 220 scaling 7.500000 only32 0
cpucycles tracesetup 1 default-perfevent precision 367 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 419 scaling 3.800000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 4130 scaling 3800.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3800000000
cpucycles implementation ppc64-mftb
cpucycles iqm 131 +26+27-4-3+4-4+4-3+4-4-3+11-3+4-4-3+4-4+4-3+4-4+4-3+4-4+4+4-3-4+4-3+4+11-3-11+4-4+4-3+4-4-3-4+12-4+4-3+4-11+4-4+4-3-4+12+4-11+4-4-3+4-4-3
cpucycles observed persecond 2535000000...5449000000 with 2048 loops 3 microseconds
cpucycles observed persecond 3042000000...5350000000 with 4096 loops 4 microseconds
cpucycles observed persecond 3244222222...4275000000 with 8192 loops 8 microseconds
cpucycles observed persecond 3587533333...4187846154 with 16384 loops 14 microseconds
cpucycles observed persecond 3714333333...3981290323 with 32768 loops 32 microseconds
cpucycles observed persecond 3693138888...3924264706 with 65536 loops 35 microseconds
cpucycles observed persecond 3795971428...3913897059 with 131072 loops 69 microseconds
cpucycles observed persecond 3821978417...3880839417 with 262144 loops 138 microseconds
cpucycles observed persecond 3821845323...3851086957 with 524288 loops 277 microseconds
cpucycles observed persecond 3835166064...3849797102 with 1048576 loops 553 microseconds
```

`cfarm45`,
AMD Athlon II X4 640,
Debian 8.11,
Linux kernel 3.16.0-11-686-pae:
```
cpucycles version 20230105
cpucycles tracesetup 0 x86-tsc precision 199 scaling 1.000000 only32 0
cpucycles tracesetup 1 x86-tscasm precision 199 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 170 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 941 scaling 3.000000 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 3200 scaling 3000.000000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3000000000
cpucycles implementation default-perfevent
cpucycles median 72 +12+0+0+0+0+0+0+0+5+0+0+0+0+0+0+0+2+0+0+0+0+0+0+0+1+0+0+0+0+0+0+0+2+0+0+0+0+0+0+0+1+0+0+0+0+0+0+0+2+0+0+0+0+0+0+0+1+0+0+0+0+0+0
cpucycles observed persecond 541500000...1812000000 with 1024 loops 3 microseconds
cpucycles observed persecond 712333333...1212250000 with 2048 loops 5 microseconds
cpucycles observed persecond 1193285714...1733600000 with 4096 loops 6 microseconds
cpucycles observed persecond 1689176470...1804562500 with 8192 loops 33 microseconds
cpucycles observed persecond 1713074626...1770600000 with 16384 loops 66 microseconds
cpucycles observed persecond 1765107692...1795140625 with 32768 loops 129 microseconds
cpucycles observed persecond 1785369649...1800603922 with 65536 loops 256 microseconds
cpucycles observed persecond 1781377862...1796288462 with 131072 loops 261 microseconds
cpucycles observed persecond 1772647398...1778247827 with 262144 loops 691 microseconds
cpucycles observed persecond 1789670493...1794149598 with 524288 loops 870 microseconds
cpucycles observed persecond 1860276211...1861561332 with 1048576 loops 3156 microseconds
```

`cfarm91`,
StarFive JH7100,
Debian trixie/sid,
Linux kernel 5.18.11-starfive:
```
cpucycles version 20250925
cpucycles tracesetup 0 riscv64-rdcycle precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 2702 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 1351 scaling 2.399988 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 2599 scaling 2399.987654 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation default-monotonic
cpucycles iqm 1476 +828-324+60+60+60+60-324+60-324+60+60-324+60+60-324+60+60-324+60-324+60+60-324+62-324+60+60-324+60+60-324+60-324+60+60-324+60+60-324+60-324+60+60-324+60+60-324+60-324+60+60+60+60+60-324+60-324+62+60+444+60-324+60+60
cpucycles observed persecond 1440250000...3968333334 with 1024 loops 7 microseconds
cpucycles observed persecond 1920181818...3072222223 with 2048 loops 10 microseconds
cpucycles observed persecond 2162473684...2733294118 with 4096 loops 18 microseconds
cpucycles observed persecond 2218777777...2530000000 with 8192 loops 35 microseconds
cpucycles observed persecond 2338014705...2490318182 with 16384 loops 67 microseconds
cpucycles observed persecond 2358559701...2440840910 with 32768 loops 133 microseconds
cpucycles observed persecond 2380905660...2419452472 with 65536 loops 264 microseconds
cpucycles observed persecond 2390815939...2410163810 with 131072 loops 526 microseconds
cpucycles observed persecond 2393901140...2403580953 with 262144 loops 1051 microseconds
cpucycles observed persecond 2397546190...2402577217 with 524288 loops 2099 microseconds
cpucycles observed persecond 2398864140...2401570967 with 1048576 loops 4327 microseconds
```

`cfarm92`,
SiFive Freedom U740,
Ubuntu 22.04.3,
Linux kernel 5.19.0-1021-generic:
```
cpucycles version 20240114
cpucycles tracesetup 0 riscv64-rdcycle precision 8 scaling 1.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 2599 scaling 2.399988 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 2599 scaling 2399.987654 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation riscv64-rdcycle
cpucycles median 8 +168+20+2+2+0+0+0+0+570+0+0+0+0+0+0+0+144+0+0+0+0+0+0+0+160+0+0+0+0+0+0+0+160+0+0+0+0+0+0+0+154+0+0+0+0+0+0+0+154+0+0+0+0+0+0+0+152+0+0+0+0+0+0
cpucycles observed persecond 571500000...2198000000 with 1024 loops 3 microseconds
cpucycles observed persecond 833600000...2094000000 with 2048 loops 4 microseconds
cpucycles observed persecond 921888888...1445142858 with 4096 loops 8 microseconds
cpucycles observed persecond 1029625000...1320642858 with 8192 loops 15 microseconds
cpucycles observed persecond 1137034482...1284481482 with 16384 loops 28 microseconds
cpucycles observed persecond 1155701754...1227454546 with 32768 loops 56 microseconds
cpucycles observed persecond 1177464285...1217163637 with 65536 loops 111 microseconds
cpucycles observed persecond 1188018099...1207858448 with 131072 loops 220 microseconds
cpucycles observed persecond 1189925170...1200519363 with 262144 loops 440 microseconds
cpucycles observed persecond 1193962457...1199117446 with 524288 loops 878 microseconds
cpucycles observed persecond 1194051324...1196780111 with 1048576 loops 1811 microseconds
```

`cfarm94`,
SiFive StarFive JH7110,
Alpine 3.22.0_alpha20250108,
Linux kernel 6.12.19-2-lts:
```
cpucycles version 20260625
cpucycles tracesetup 0 riscv64-perfrdcycle precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 riscv64-rdcycle precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 0 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 1750 scaling 1.500000 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 2510 scaling 1500.000000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1500000000
cpucycles implementation default-monotonic
cpucycles iqm 1125 +0+0+0-375+0+0+0+0+0+0-375+0+0+0-375+0+0+0+0-375+0+0+0-375+0+0+0-375+0+0+0+0-375+0+0+0-375+0+0+0+0-375+0+0+1-375+0+0+0-375+0+0+0+0-375+0+0+0-375+0+0+0-375+0
cpucycles observed persecond 600000000...2500000000 with 1024 loops 4 microseconds
cpucycles observed persecond 875000000...2343750000 with 2048 loops 5 microseconds
cpucycles observed persecond 1083444444...1982285715 with 4096 loops 8 microseconds
cpucycles observed persecond 1258928571...1812500000 with 8192 loops 13 microseconds
cpucycles observed persecond 1312576923...1609458334 with 16384 loops 25 microseconds
cpucycles observed persecond 1412255319...1583355556 with 32768 loops 46 microseconds
cpucycles observed persecond 1451762376...1526545455 with 65536 loops 100 microseconds
cpucycles observed persecond 1487322033...1527891429 with 131072 loops 176 microseconds
cpucycles observed persecond 1490735537...1511465374 with 262144 loops 362 microseconds
cpucycles observed persecond 1493099573...1503243938 with 524288 loops 702 microseconds
cpucycles observed persecond 1497124823...1502417373 with 1048576 loops 1417 microseconds
```

`cfarm95`,
SpacemiT X60,
Debian trixie/sid,
Linux kernel 6.6.36-cfarm #1:
```
cpucycles version 20260625
cpucycles tracesetup 0 riscv64-perfrdcycle precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 riscv64-rdcycle precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-perfevent precision 1794 scaling 1.000000 only32 0
cpucycles tracesetup 3 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-monotonic precision 1111 scaling 1.228800 only32 0
cpucycles tracesetup 5 default-gettimeofday precision 2248 scaling 1228.800000 only32 0
cpucycles tracesetup 6 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 1228800000
cpucycles implementation default-monotonic
cpucycles iqm 135 -33+19-32+172+19-33+18-33+19+19-32+888+19-33+19-32+18-33+19+224-33+18-33+19-32+19-33+377-33+19-33+18-31+18-33+428-32+19-33+18-33+19-32+633+19-33+18-33+19-32+19+376+18-31+18-33+19-33+18+19+19-32+19-33
cpucycles observed persecond 468000000...6052000000 with 1024 loops 6 microseconds
cpucycles observed persecond 962800000...2696333334 with 2048 loops 4 microseconds
cpucycles observed persecond 1043250000...1476166667 with 4096 loops 7 microseconds
cpucycles observed persecond 1102600000...1311461539 with 8192 loops 14 microseconds
cpucycles observed persecond 1175785714...1285923077 with 16384 loops 27 microseconds
cpucycles observed persecond 1198090909...1253943397 with 32768 loops 54 microseconds
cpucycles observed persecond 1214266055...1246056075 with 65536 loops 108 microseconds
cpucycles observed persecond 1217488479...1237860466 with 131072 loops 216 microseconds
cpucycles observed persecond 1224293023...1233245328 with 262144 loops 429 microseconds
cpucycles observed persecond 1225245348...1230967366 with 524288 loops 859 microseconds
cpucycles observed persecond 1227769590...1230615926 with 1048576 loops 1709 microseconds
```

`cfarm103`,
Apple M1 (Icestorm-M1 + Firestorm-M1),
Debian trixie/sid,
Linux kernel 6.5.0-asahi-00780-g62806c2c6f29:
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1086 scaling 86.000000 only32 0
cpucycles tracesetup 3 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1094 scaling 2.064000 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 3084 scaling 2064.000000 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2064000000
cpucycles implementation arm64-vct
cpucycles iqm 0 +0+0+0+86+0+0+0+0+258+0+0+0+0+0+0+0+86+0+0+0+0+0+0+0+0+86+0+0+0+0+0+0+0+0+86+0+0+0+0+0+0+0+0+86+0+0+0+0+0+0+0+0+86+0+0+0+0+0+0+0+86+0+0+0
cpucycles observed persecond 1311500000...2709000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1748666666...2666000000 with 8192 loops 5 microseconds
cpucycles observed persecond 1834666666...2383428572 with 16384 loops 8 microseconds
cpucycles observed persecond 2053250000...2358857143 with 32768 loops 15 microseconds
cpucycles observed persecond 1988424242...2122258065 with 65536 loops 32 microseconds
cpucycles observed persecond 2017692307...2084476191 with 131072 loops 64 microseconds
cpucycles observed persecond 2050129032...2078052288 with 262144 loops 154 microseconds
cpucycles observed persecond 2056833333...2071818182 with 524288 loops 287 microseconds
cpucycles observed persecond 2060477815...2067828768 with 1048576 loops 585 microseconds
```

`cfarm104`,
Apple M1 (Icestorm-M1 + Firestorm-M1),
MacOSX 12.6 21.6.0:
```
cpucycles version 20251226
cpucycles tracesetup 0 arm64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 200 scaling 100.000000 only32 0
cpucycles tracesetup 3 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 4 default-mach precision 200 scaling 100.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 2699 scaling 2.399988 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 2699 scaling 2399.987654 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation arm64-vct
cpucycles iqm 0 +4500+0+0+0+0+0+100+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+100+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+100+0+0+0+0+0+0+0
cpucycles observed persecond 1725000000...3450000000 with 8192 loops 3 microseconds
cpucycles observed persecond 1885714285...2640000000 with 16384 loops 6 microseconds
cpucycles observed persecond 2300000000...2533333334 with 32768 loops 22 microseconds
cpucycles observed persecond 2339285714...2523076924 with 65536 loops 27 microseconds
cpucycles observed persecond 2348529411...2421212122 with 131072 loops 67 microseconds
cpucycles observed persecond 2377678571...2420909091 with 262144 loops 111 microseconds
cpucycles observed persecond 2387214611...2410138249 with 524288 loops 218 microseconds
cpucycles observed persecond 2397814207...2411813187 with 1048576 loops 365 microseconds
```

`cfarm110` (`gcc1-power7`),
IBM POWER7,
CentOS 7.9 AltArch,
Linux kernel 3.10.0-1160.105.1.el7.ppc64:
```
cpucycles version 20260105
cpucycles tracesetup 0 ppc64-mftb precision 212 scaling 7.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 236 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 465 scaling 3.550000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 3880 scaling 3550.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3550000000
cpucycles implementation ppc64-mftb
cpucycles iqm 167 +36-34+36-34-13+15-20+22-27+29-41+43-48+36-41+36-34+22-20+15-13+8-6-6+8-13+15-20+22-34+36+15-13+8-6+1+1-13+15-20+22-34+36-48+36-41+8-6+36-34+22-20+15-13+8-6-6+8-20+22-27+29-41+36
cpucycles observed persecond 3159333333...4905250000 with 4096 loops 5 microseconds
cpucycles observed persecond 3386727272...4216333334 with 8192 loops 10 microseconds
cpucycles observed persecond 3348333333...3739842106 with 16384 loops 20 microseconds
cpucycles observed persecond 3534829268...3733153847 with 32768 loops 40 microseconds
cpucycles observed persecond 3473016129...3599166667 with 65536 loops 61 microseconds
cpucycles observed persecond 3541000000...3609900000 with 131072 loops 111 microseconds
cpucycles observed persecond 3571450450...3606177273 with 262144 loops 221 microseconds
cpucycles observed persecond 3578959367...3596333334 with 524288 loops 442 microseconds
cpucycles observed persecond 3578793453...3587492082 with 1048576 loops 885 microseconds
```

`cfarm112` (`gcc2-power8`),
IBM POWER8E,
CentOS 7.9 AltArch,
Linux kernel 3.10.0-1127.13.1.el7.ppc64le:
```
cpucycles version 20260625
cpucycles tracesetup 0 ppc64-mftb precision 1159 scaling 7.250000 only32 0
cpucycles tracesetup 1 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 1342 scaling 3.690000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 4710 scaling 3690.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3690000000
cpucycles implementation ppc64-mftb
cpucycles iqm 194 +9-6-5+1-5+9-6+9-12+1+9+2+2+2+9-13+2+16-13-5+1+2+2-6+17-6-5+9-6-5+1+2+9+2+9-6-5+16-5+1-5+9+2+1+2-13+9+9-5-6+9+2-5+1-5+9-6+17-6-13+9+9-5+2
cpucycles observed persecond 2802000000...9240500000 with 1024 loops 3 microseconds
cpucycles observed persecond 2977714285...4412200000 with 2048 loops 6 microseconds
cpucycles observed persecond 3208083333...3927300000 with 4096 loops 11 microseconds
cpucycles observed persecond 3486608695...3855952381 with 8192 loops 22 microseconds
cpucycles observed persecond 3606500000...3832027778 with 16384 loops 37 microseconds
cpucycles observed persecond 3644816666...3781258621 with 32768 loops 59 microseconds
cpucycles observed persecond 3672697478...3740752137 with 65536 loops 118 microseconds
cpucycles observed persecond 3687067510...3721102128 with 131072 loops 236 microseconds
cpucycles observed persecond 3702415254...3719497873 with 262144 loops 471 microseconds
cpucycles observed persecond 3705988335...3714527099 with 524288 loops 942 microseconds
cpucycles observed persecond 3708108821...3712812798 with 1048576 loops 1892 microseconds
```

`cfarm120`,
IBM POWER10,
AlmaLinux 9.7,
Linux kernel 5.14.0-427.31.1.el9_4.ppc64le:
```
cpucycles version 20260625
cpucycles tracesetup 0 ppc64-mftb precision 1023 scaling 5.750000 only32 0
cpucycles tracesetup 1 default-perfevent precision 139 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 1046 scaling 2.950000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 3980 scaling 2950.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2950000000
cpucycles implementation default-perfevent
cpucycles iqm 136 +45+44+0+1+2+1+0+0+0+0+0+5+2+1+3+0+0+1+1+14+1+0+1+0+0+2+0+0+0+2+0+0+0+14+0+0+1+0+1+2+1+0+0+14+0+1+0+14+1+0+1+0+0+1+3+0+1+2+0+0+0+0+0+1
cpucycles observed persecond 1811000000...3992000000 with 2048 loops 3 microseconds
cpucycles observed persecond 2836600000...4901333334 with 4096 loops 4 microseconds
cpucycles observed persecond 2868444444...3754428572 with 8192 loops 8 microseconds
cpucycles observed persecond 3355777777...3810062500 with 16384 loops 17 microseconds
cpucycles observed persecond 3409464285...3692000000 with 32768 loops 27 microseconds
cpucycles observed persecond 3736972222...3969617648 with 65536 loops 35 microseconds
cpucycles observed persecond 3785605633...3901246377 with 131072 loops 70 microseconds
cpucycles observed persecond 3830314285...3888688406 with 262144 loops 139 microseconds
cpucycles observed persecond 3873345323...3902855073 with 524288 loops 277 microseconds
cpucycles observed persecond 3884176795...3899293901 with 1048576 loops 542 microseconds
```

`cfarm185`,
Ampere eMAG 8180,
AlmaLinux 8.10,
Linux kernel 4.18.0-553.el8_10.aarch64:
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1075 scaling 75.000000 only32 0
cpucycles tracesetup 3 default-perfevent precision 205 scaling 1.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1167 scaling 3.000000 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 4030 scaling 3000.000000 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3000000000
cpucycles implementation default-perfevent
cpucycles iqm 264 -25-16-17+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 890000000...4520500000 with 1024 loops 3 microseconds
cpucycles observed persecond 1082000000...1897250000 with 2048 loops 5 microseconds
cpucycles observed persecond 1811000000...2751800000 with 4096 loops 6 microseconds
cpucycles observed persecond 2257272727...2859777778 with 8192 loops 10 microseconds
cpucycles observed persecond 2607736842...2970588236 with 16384 loops 18 microseconds
cpucycles observed persecond 2901970588...3114531250 with 32768 loops 33 microseconds
cpucycles observed persecond 3028430769...3138857143 with 65536 loops 64 microseconds
cpucycles observed persecond 3074007812...3128968254 with 131072 loops 127 microseconds
cpucycles observed persecond 3121726190...3149980000 with 262144 loops 251 microseconds
cpucycles observed persecond 3127475149...3141540919 with 524288 loops 502 microseconds
cpucycles observed persecond 3142828171...3149926927 with 1048576 loops 1000 microseconds
```

`cfarm202`,
UltraSparc T5,
Debian forky/sid,
Linux kernel 6.17.0-rc5+:
```
cpucycles version 20260625
cpucycles tracesetup 0 sparc64-rdtick precision 65 scaling 1.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 448 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 1235 scaling 3.599910 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 4629 scaling 3599.910000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 3599910000
cpucycles implementation sparc64-rdtick
cpucycles iqm 73 +24+0+24+24+24+0+0+24+0+0+1+0+0+0+0+0+0+0+0+0+0+1+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 2888000000...4979666667 with 4096 loops 4 microseconds
cpucycles observed persecond 3197777777...4196285715 with 8192 loops 8 microseconds
cpucycles observed persecond 3380235294...3863000000 with 16384 loops 16 microseconds
cpucycles observed persecond 3479030303...3719000000 with 32768 loops 32 microseconds
cpucycles observed persecond 3530707692...3650428572 with 65536 loops 64 microseconds
cpucycles observed persecond 3584937500...3645658731 with 131072 loops 127 microseconds
cpucycles observed persecond 3584468750...3614586615 with 262144 loops 255 microseconds
cpucycles observed persecond 3591248532...3606304519 with 524288 loops 510 microseconds
cpucycles observed persecond 3594648383...3602175663 with 1048576 loops 1020 microseconds
```

`cfarm216`,
VM on SPARC-M8,
Solaris 11.4.88.207.1:
```
cpucycles version 20260105
cpucycles tracesetup 0 sparc64-rdtick precision 63 scaling 1.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 614 scaling 5.067000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 5387 scaling 5067.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 5067000000
cpucycles implementation sparc64-rdtick
cpucycles iqm 69 +0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+57376+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 2865000000...5940000000 with 4096 loops 3 microseconds
cpucycles observed persecond 4597800000...7831000000 with 8192 loops 4 microseconds
cpucycles observed persecond 4558600000...5760750000 with 16384 loops 9 microseconds
cpucycles observed persecond 4765578947...5356647059 with 32768 loops 18 microseconds
cpucycles observed persecond 4794418604...5041268293 with 65536 loops 42 microseconds
cpucycles observed persecond 4968216216...5114986112 with 131072 loops 73 microseconds
cpucycles observed persecond 4965598484...5044946154 with 262144 loops 131 microseconds
cpucycles observed persecond 4984102661...5023766284 with 524288 loops 262 microseconds
cpucycles observed persecond 4985677238...5005539326 with 1048576 loops 535 microseconds
```

`cfarm230`,
Cavium Octeon III V0.2,
Debian 10.13,
Linux kernel 4.9.79-UBNT_E300:
```
cpucycles version 20260625
cpucycles tracesetup 0 mips64-cc precision 323 scaling 1.000000 only32 1
cpucycles tracesetup 1 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 2243 scaling 2.399988 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 3419 scaling 2399.987654 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation mips64-cc
cpucycles iqm 405 +30+17+6+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0
cpucycles observed persecond 868444444...1931285715 with 1024 loops 8 microseconds
cpucycles observed persecond 923437500...1160571429 with 2048 loops 15 microseconds
cpucycles observed persecond 969766666...1089250000 with 4096 loops 29 microseconds
cpucycles observed persecond 979067796...1037964913 with 8192 loops 58 microseconds
cpucycles observed persecond 983837606...1013113044 with 16384 loops 116 microseconds
cpucycles observed persecond 994792207...1009589520 with 32768 loops 230 microseconds
cpucycles observed persecond 996036876...1003851852 with 65536 loops 460 microseconds
cpucycles observed persecond 998830250...1002534352 with 131072 loops 918 microseconds
cpucycles observed persecond 998601196...1000898149 with 262144 loops 1837 microseconds
cpucycles observed persecond 999575323...1000699892 with 524288 loops 3703 microseconds
cpucycles observed persecond 999797523...1000292365 with 1048576 loops 7348 microseconds
```

`cfarm400`,
Loongson-3C5000L-LL,
Debian trixie/sid,
Linux kernel 6.1.0-rc7+:
```
cpucycles version 20260105
cpucycles tracesetup 0 loong64-rdtime precision 124 scaling 24.000000 only32 0
cpucycles tracesetup 1 default-perfevent precision 180 scaling 1.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 416 scaling 2.399988 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 2729 scaling 2399.987654 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation loong64-rdtime
cpucycles iqm 48 +24+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+24+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+24+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+24+0+0+0+0+0+0+0
cpucycles observed persecond 1982400000...3408000000 with 4096 loops 4 microseconds
cpucycles observed persecond 1975200000...2514000000 with 8192 loops 9 microseconds
cpucycles observed persecond 2190666666...2482500000 with 16384 loops 17 microseconds
cpucycles observed persecond 2315294117...2469000000 with 32768 loops 33 microseconds
cpucycles observed persecond 2348776119...2425846154 with 65536 loops 66 microseconds
cpucycles observed persecond 2383636363...2422892308 with 131072 loops 131 microseconds
cpucycles observed persecond 2392425855...2411954023 with 262144 loops 262 microseconds
cpucycles observed persecond 2392334600...2402061069 with 524288 loops 525 microseconds
cpucycles observed persecond 2396822857...2401694657 with 1048576 loops 1049 microseconds
```

`cfarm425`,
Ampere Altra Max M128-30,
Debian forky,
Linux kernel 6.16.12+deb14+1-arm64:
```
cpucycles version 20260625
cpucycles tracesetup 0 arm64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 arm64-pmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 arm64-vct precision 1096 scaling 96.000000 only32 0
cpucycles tracesetup 3 default-perfevent precision 91 scaling 1.000000 only32 0
cpucycles tracesetup 4 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 5 default-monotonic precision 1116 scaling 2.399988 only32 0
cpucycles tracesetup 6 default-gettimeofday precision 3429 scaling 2399.987654 only32 0
cpucycles tracesetup 7 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation default-perfevent
cpucycles iqm 82 +29-1-1-1-1-1+386-1-1-1+179-1-1+130+161-1+235-1-1-1-1+153+170-1-1+149-1-1-1+258+154-1-1-1-1-1-1-1-1-1-1+196+194-1+172-1-1-1-1-1-1-1+158-1-1+183-1-1-1-1-1-1-1+235
cpucycles observed persecond 547750000...1263500000 with 2048 loops 3 microseconds
cpucycles observed persecond 852200000...1532333334 with 4096 loops 4 microseconds
cpucycles observed persecond 1412166666...2203250000 with 8192 loops 5 microseconds
cpucycles observed persecond 1822600000...2320375000 with 16384 loops 9 microseconds
cpucycles observed persecond 2064187500...2383285715 with 32768 loops 15 microseconds
cpucycles observed persecond 2266612903...2434482759 with 65536 loops 30 microseconds
cpucycles observed persecond 2322280701...2412890910 with 131072 loops 56 microseconds
cpucycles observed persecond 2391572727...2438981482 with 262144 loops 109 microseconds
cpucycles observed persecond 2413564220...2437481482 with 524288 loops 217 microseconds
cpucycles observed persecond 2432299539...2444340278 with 1048576 loops 433 microseconds
```

`cfarm430`,
AMD EPYC 7773X,
FreeBSD 16.0 VM:
```
cpucycles version 20260625
cpucycles tracesetup 0 amd64-perfpmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 1 amd64-perfpmc precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 amd64-pmcff precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 amd64-tsc precision 1022 scaling 1.000000 only32 0
cpucycles tracesetup 4 amd64-tscasm precision 1032 scaling 1.000000 only32 0
cpucycles tracesetup 5 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 6 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 7 default-monotonic precision 1140 scaling 2.399988 only32 0
cpucycles tracesetup 8 default-gettimeofday precision 3429 scaling 2399.987654 only32 0
cpucycles tracesetup 9 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 2399987654
cpucycles implementation amd64-tsc
cpucycles iqm 26 +40-4-4+18-4-4+150+18-4+18-4-4-4+18-4-4+18-4+18-4-4+18-4-4-4+18-4-4+18-4-4+18-4-4+18-4-4-4+18-4-4+18-4+18-4-4+18-4-4+18-4-4-4+18-4-4+18-4+18-4-4+18-4-4
cpucycles observed persecond 1578500000...3322000000 with 4096 loops 3 microseconds
cpucycles observed persecond 1782000000...2569600000 with 8192 loops 6 microseconds
cpucycles observed persecond 2053333333...2503600000 with 16384 loops 11 microseconds
cpucycles observed persecond 2053040000...2245913044 with 32768 loops 24 microseconds
cpucycles observed persecond 2189423076...2283600000 with 65536 loops 51 microseconds
cpucycles observed persecond 2179614678...2224467290 with 131072 loops 108 microseconds
cpucycles observed persecond 2189704433...2216746269 with 262144 loops 202 microseconds
cpucycles observed persecond 2195588235...2211413979 with 524288 loops 373 microseconds
cpucycles observed persecond 2196466750...2203708177 with 1048576 loops 796 microseconds
```

`z15`,
IBM z15:
```
cpucycles version 20230106
cpucycles tracesetup 0 s390x-stckf precision 250 scaling 1.269531 only32 0
cpucycles tracesetup 1 default-perfevent precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 2 default-mach precision 0 scaling 0.000000 only32 0
cpucycles tracesetup 3 default-monotonic precision 272 scaling 5.200000 only32 0
cpucycles tracesetup 4 default-gettimeofday precision 5400 scaling 5200.000000 only32 0
cpucycles tracesetup 5 default-zero precision 0 scaling 0.000000 only32 0
cpucycles persecond 5200000000
cpucycles implementation s390x-stckf
cpucycles median 48 +87+8+0-2+0+0+38-2+0+1-3+1+28+0+3-3+1+0+28+0-2+3+0-2+36+0+0+0+1+0+28+0-2+0+3-2+35+1+0-2+0+3+28+0-2+0+0-2+3+25+3+0-2+0+1+35+1+0+0-2+0+28+0
cpucycles observed persecond 4948941176...5627733334 with 8192 loops 16 microseconds
cpucycles observed persecond 4104125000...5515666667 with 16384 loops 7 microseconds
cpucycles observed persecond 5047076923...5987818182 with 32768 loops 12 microseconds
cpucycles observed persecond 5044846153...5475708334 with 65536 loops 25 microseconds
cpucycles observed persecond 5141313725...5357428572 with 131072 loops 50 microseconds
cpucycles observed persecond 5150892156...5257250000 with 262144 loops 101 microseconds
cpucycles observed persecond 5183421568...5236549505 with 524288 loops 203 microseconds
cpucycles observed persecond 5190282555...5216582717 with 1048576 loops 406 microseconds
```
