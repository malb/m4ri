To download and unpack the latest version of libcpucycles:

    wget -m https://cpucycles.cr.yp.to/libcpucycles-latest-version.txt
    version=$(cat cpucycles.cr.yp.to/libcpucycles-latest-version.txt)
    wget -m https://cpucycles.cr.yp.to/libcpucycles-$version.tar.gz
    tar -xzf cpucycles.cr.yp.to/libcpucycles-$version.tar.gz
    cd libcpucycles-$version

Then [install](install.html).

### Archives and changelog (reverse chronological)

[`libcpucycles-20260625.tar.gz`](libcpucycles-20260625.tar.gz) [browse](libcpucycles-20260625.html)

In `amd64-perfpmcff` and `default-perfevent`,
to compensate for the Linux kernel
[not doing](https://cr.yp.to/2026/20260624-rdpmc.c)
what the `perf_event_open` manual page promises,
try opening a second fd using PMU type 10 for Intel E-cores.
This should prevent fallbacks to `tsc` for programs pinned to E-cores,
and should let cycle counting continue for programs that migrate from P-cores to E-cores.

Clear `/proc/sys/kernel/nmi_watchdog` in `cpucycles-open`.
Otherwise the Linux NMI watchdog forces a wraparound after some number of seconds
in RDPMC cycle counts from `amd64-perfpmcff`,
corrupting multi-second cycle counts on Intel machines.
AMD machines are unaffected.

Increase bump to 1000 for all off-core counters and fixed-frequency counters.

---

[`libcpucycles-20260105.tar.gz`](libcpucycles-20260105.tar.gz) [browse](libcpucycles-20260105.html)

Bigger split of RDPMC-based counters
to handle not just `perf` sometimes incorrectly returning index 0 on Intel
but also AMD not supporting fixed-function counter `0x40000001`.
New `amd64-perfpmcff` tries fixed-function counter `0x40000001` within `perf`;
`amd64-perfpmc` now goes back to trying index returned by `perf`;
`amd64-pmcff` (renamed from `amd64-pmc`) tries fixed-function counter `0x40000001` outside `perf`.

Add `orderbump` mechanism,
for example to prefer `tsc` more systematically over `tscasm`.

---

[`libcpucycles-20251226.tar.gz`](libcpucycles-20251226.tar.gz) [browse](libcpucycles-20251226.html)

Split `amd64-pmc`, `arm64-pmc`, `riscv32-rdcycle`, and `riscv64-rdcycle`
into `perf` and non-`perf` versions,
with the `perf` versions tried first.
Hopefully the `perf` versions will also eliminate the need to install kernel modules
for `amd64-pmc` and `arm64-pmc`.

For `amd64-pmc`,
always use fixed-function cycle counter `(1<<30)|1`.
This works around a bug in, e.g., kernel 6.1.0 `perf` incorrectly returning index 0 for E-cores.

For all `perf`-based counters,
use `PERF_FLAG_FD_CLOEXEC` if it is available.

For `*tsc*` counters,
return invariant TSC frequency if that is clear from CPUID information.
Thanks to Tee-Kiah Chia for the suggestion and a prototype patch.

Add internal `ticks_close()`,
most importantly for closing `perf` descriptors.
(An alternative would be to unify the `perf` handling
inside a single `perf_event_open` in `wrapper` conditioned on Linux,
but this raises the question of how to handle `perf` platform variations
such as `arm64-pmc` needing `3` in `config1`.)

Increase penalty for `EXTEND32` counters.
(It might be better to directly benchmark the extension.)

Increase penalty for fixed-frequency counters,
including RDTSC when it is known to be invariant.

For `wrapper`,
use `sigjmp_buf` instead of `jmp_buf`.

Document requirements on
`/sys/devices/cpu*/rdpmc` and
`/proc/sys/kernel/perf_user_access`.
Add `cpucycles-open` abstraction.

Update `configure`:
recognize a few more x86 synonyms;
`--no-` instead of `--no`;
support `configure.log`.

Support `$DESTDIR` in `install`.
Thanks to Robert Clausecker for the suggestion.

---

[`libcpucycles-20250925.tar.gz`](libcpucycles-20250925.tar.gz) [browse](libcpucycles-20250925.html)

For `perfevent`,
try all combinations of disabling exclude-kernel and disabling exclude-hv,
to handle platforms where cycle counters cannot handle one or both of those.

For `perfevent`,
do `DISABLE` and `ENABLE` around counter read (suggested by Jim Apple)
to (hopefully) handle Graviton 3.
But skip this if a setup test seems to work without it.

Recognize `sun4v` as `sparc64`.

Add support for `loong64`.

Eliminate various compiler warnings.

For `cpucycles-info`,
print interquartile mean (iqm) rather than median.

Say "sometimes" in documentation for enabling cycle counters in ARM via kernel modules.

Change HTML style, in particular for better usability on phones.

---

[`libcpucycles-20240318.tar.gz`](libcpucycles-20240318.tar.gz) [browse](libcpucycles-20240318.html)

Port to MacOS X:
handle missing `-lrt`, and handle differences in shared-library naming.

Include `cpucycles-info` man page.

---

[`libcpucycles-20240114.tar.gz`](libcpucycles-20240114.tar.gz) [browse](libcpucycles-20240114.html)

Add `arm32-1176` counter.

Allow slop 0.2 rather than 0.1 for `FINDMULTIPLIER`.

Improve platform detection.

Port to FreeBSD.

Use blue boldface during compilation for "skipping option that did not compile".

`doc/install.md`: headings; note manual pages.

Add `doc/license.md`.

Update HTML style for better tt visibility and copy-paste.

---

[`libcpucycles-20230115.tar.gz`](libcpucycles-20230115.tar.gz) [browse](libcpucycles-20230115.html)

Update actual `cpucycles_version` behavior to match documentation.

---

[`libcpucycles-20230110.tar.gz`](libcpucycles-20230110.tar.gz) [browse](libcpucycles-20230110.html)

`doc/api.md`: Document `cpucycles_version()`.

Add `s390x-stckf` counter.

`cpucycles/default-perfevent.c`: Read into `int64_t` instead of `long long`.
Add comment explaining issues with `PERF_FORMAT_TOTAL_TIME_RUNNING`.

`configure`: Improve `uname` handling.

`doc/api.md`: Update description of default frequency.

---

[`libcpucycles-20230105.tar.gz`](libcpucycles-20230105.tar.gz) [browse](libcpucycles-20230105.html)

Initial release.
