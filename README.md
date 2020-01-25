M4RI is a library for fast arithmetic with dense matrices over F2. The name M4RI comes from the first implemented algorithm: The “Method of the Four Russians” inversion algorithm published by Gregory Bard. This algorithm in turn is named after the “Method of the Four Russians” multiplication algorithm which is probably better referred to as Kronrod's method. M4RI is available under the General Public License Version 2 or later (GPLv2+).

# Main Features #

* basic arithmetic with dense matrices over F2 (addition, equality testing, stacking, augmenting, sub-matrices, randomisation);

* asymptotically fast $O(n^{log_2 7})$ matrix multiplication via the Method of the Four Russians (M4RM) & Strassen-Winograd algorithm;

* asymptotically fast $O(n^{log_2 7})$ PLE factorisation (Gaussian elimination, system solving, …);

* fast row echelon form computation and matrix inversion via the Method of the Four Russians (M4RI, $O(n^{3/log n})$);

* asymptotically fast Triangular System solving with Matrices (upper left, lower left, upper right, lower right),

* support for the x86/x86_64 SSE2 instruction set where available;

* preliminary support for parallelisation on shared memory systems via OpenMP;

* and support for Linux, Solaris, and OS X (GCC).

See [Further Reading](https://bitbucket.org/malb/m4ri/wiki/Further%20Reading) for implemented algorithms.

# Performance

See [Performance](https://bitbucket.org/malb/m4ri/wiki/Performance).

# OpenMP Support #

OpenMP support for parallel multiplication and elimination is enabled with the

    --enable-openmp

configure switch.

# Install #

If you downloaded M4RI by cloning the mainline tree at

https://bitbucket.org/malb/m4ri

you need to first run the following command:

    autoreconf --install

Then do the usual

    ./configure
    make
    make check

For details see the instructions in the file `INSTALL`.

# Documentation #

To build the reference manual, ensure that you have Doxygen installed. The HTML version of the reference manual can be built as follows:

    cd src/
    doxygen

The built documentation is contained under the doc subdirectory of m4ri/. Once the HTML version is built, you can build the PDF version as follows:

    cd doc/latex/
    make

The documentation is also available [here](http://malb.bitbucket.io/m4ri/).

# Contributors

At least the following people have contributed to the M4RI library.

* **[Tim Abbott](http://web.mit.edu/tabbott/www/)**: Debian-isation & advice on correct libtool versioning;
* **[Martin Albrecht](http://martinralbrecht.wordpress.com)**: maintainer, release manager, peformance tuning (M4RM, M4RI, Strassen, PLE), initial M4RM implementation, parallelisation, PLE factorisation (MMPF algorithm);
* **[Gregory Bard](http://www.math.umd.edu/~bardg/)**: initial author, M4RI algorithm and initial implementation;
* **[Marco Bodrato](http://bodrato.it/)**: new [Strassen-like sequence](http://bodrato.it/software/strassen.html) for matrix multiplication and squaring which improves performance for squaring;
* **[Michael Brickenstein](http://www.mfo.de/organisation/institute/brickenstein/)**: [PolyBoRi](http://polybori.sourceforge.net) author, standard conformity contributions for ANSIC, test data, discussion/suggestion of performance improvements, fast vector-matrix products;
* **[Alexander Dreyer](http://www.itwm.fhg.de/en/as__asemployees__dreyer/dreyer/)**: [PolyBoRi](http://polybori.sourceforge.net) author, standard conformity contributions for ANSIC;
* **[Jean-Guillaume Dumas](http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/)**: linear system resolution;
* **[William Hart](http://www.warwick.ac.uk/~masfaw/)**: many performance improvements for matrix multiplication and in general;
* **[David Harvey](http://cims.nyu.edu/~harvey/)**: parallel parity function used in classical multiplication;
* **Jerry James**: bug fixes, dealing with compiler warnings, Fedora Linux packaging;
* **David Kirkby**: portability issues (Solaris, HP Unix);
* **[Clément Pernet](http://www.math.washington.edu/~pernet/)**: PLE factorisation, triangular
    system solving (TRSM);
* **Wael Said**: test cases, feedback;
* **Carlo Wood**: bit-level optimisation (transpose, column swaps), refactoring, benchmark(et)ing framework, test code, build system clean-up;

We are grateful to **[William Stein](http://modular.math.washington.edu/)** for providing our hosting and general infrastructure in the past.

# Citing M4RI

If you use our libraries in a non-trivial part of your research please consider citing them as follows:

	@manual{M4RI,
	    key          = "M4RI",
	    author       = "Martin Albrecht and Gregory Bard",
	    organization = "The M4RI~Team",
	    title        = "{The M4RI Library -- Version **version**}",
	    year         = **year**,
	    url          = "\url{http://m4ri.sagemath.org}",
	}

and cite the appropriate publications mentioned in [Further Reading](https://bitbucket.org/malb/m4ri/wiki/Further%20Reading).


# Contact

Please contact our [mailinglist](http://groups.google.com/group/m4ri-devel) if there are bugs, questions, comments.

# History

* **2020/01/25** A new version of M4RI is available with a few bugfixes. It is available at https://bitbucket.org/malb/m4ri/downloads.

* **2020/01/15** A new version of M4RI is available with a few small build system tweaks. It is available at https://bitbucket.org/malb/m4ri/downloads.

* **2015/04/17** Our hosting for http://m4ri.sagemath.org at University of Washington. is discontinued and we’re moving everything over to https://bitbucket.org/malb/m4ri. A copy of the old website (except for large files) is available at http://malb.bitbucket.io/m4ri-e-website-2008-2015/.

* **2014/09/14** A new version of M4RI and M4RIE is available for [download](https://bitbucket.org/malb/m4ri/downloads). The biggest change is that `A->offset` was dropped. Also, various small (multicore) performance improvements were implemented. The update for M4RIE is to maintain compatibility with M4RI. A few improvements were implemented for the mzd_poly module as well.

* **2013/04/16** A new version of M4RI is available for [download](https://bitbucket.org/malb/m4ri/downloads). A detailed changlog is available [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20130416) for M4RI.

* **2012/12/21** A new version of M4RI is available for [download](https://bitbucket.org/malb/m4ri/downloads). A detailed changlog is available [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20121224) for M4RI. See also this [blog post](https://martinralbrecht.wordpress.com/2012/12/21/m4ri-20121224/) for details.

* **2012/06/13** New versions of both M4RI and M4RIE are available for [download](https://bitbucket.org/malb/m4ri/downloads). A detailed changlog are available [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20120613) for M4RI.

* **2012/04/13** New versions of both M4RI and M4RIE are available for [download](https://bitbucket.org/malb/m4ri/downloads). Detailed changlogs are available [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20120415) for M4RI and [here](https://bitbucket.org/malb/m4rie/wiki/M4RIE-20120415) for M4RIE.

* **2011/12/04** New versions of both M4RI and M4RIE are available for [download](https://bitbucket.org/malb/m4ri/downloads). The highlight of this version for M4RI is support for reading and writing 1-bit PNG images. The highlight of this release of M4RIE is much improved performance for $4 < e \leq 8$. Detailed changlogs are available [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20111203) for M4RI and [here](https://bitbucket.org/malb/m4rie/wiki/M4RIE-20111203) for M4RIE.

* **2011/11/30** A [technical report](http://arxiv.org/abs/1111.6900) by Martin R. Albrecht is available describing the M4RIE library. In particular, Newton-John tables are introduced and our implementation of Karatsuba based matrix-matrix multiplication is described:

  > **The M4RIE library for dense linear algebra over small fields with even characteristic**
  >  
  > *Abstract:* In this work, we present the M4RIE library which implements efficient algorithms for
  > linear algebra with dense matrices over GF(2^e) for 2 ≤ e ≤ 10. As the name of the library
  > indicates, it makes heavy use of the M4RI library both directly (i.e., by calling it) and
  > indirectly (i.e., by using its concepts). We provide an open-source GPLv2+ C library for
  > efficient linear algebra over GF(2^e) for e small. In this library we implemented an idea due to
  > Bradshaw and Boothby which reduces matrix multiplication over GF(p^k) to a series of matrix
  > multiplications over GF(p). Furthermore, we propose a caching technique - Newton-John tables -
  > to avoid finite field multiplications which is inspired by Kronrod's method ("M4RM") for matrix
  > multiplication over GF(2). Using these two techniques we provide asymptotically fast triangular
  > solving with matrices (TRSM) and PLE-based Gaussian elimination. As a result, we are able to
  > significantly improve upon the state of the art in dense linear algebra over $F(2^e) with 2 ≤ e
  > ≤ 10.

* **2011/11/29** A [technical report](http://arxiv.org/abs/1111.6549) by Martin R. Albrecht, Gregory Bard and Clément Pernet is available describing the Gaussian elimination machinery (PLE decomposition) in the M4RI library:

  > **Efficient Dense Gaussian Elimination over the Finite Field with Two Elements.**
  >  
  > *Abstract:* In this work we describe an efficient implementation of a hierarchy of algorithms
  > for Gaussian elimination upon dense matrices over the field with two elements. We discuss both
  > well-known and new algorithms as well as our implementations in the M4RI library, which has been
  > adopted into Sage. The focus of our discussion is a block iterative algorithm for PLE
  > decomposition which is inspired by the M4RI algorithm. The implementation presented in this work
  > provides considerable performance gains in practice when compared to the previously fastest
  > implementation. We provide performance figures on x86_64 CPUs to demonstrate the alacrity of our
  > approach.

* **2011/10/10** A new release of M4RI is available for [download](https://bitbucket.org/malb/m4ri/downloads/m4ri-20111004.tar.gz). See the [release notes](https://bitbucket.org/malb/m4ri/wiki/M4RI-20111004) for the list of changes. Also, a new release of M4RIE is also available for [download](https://bitbucket.org/malb/m4ri/downloads/m4rie-20111004.tar.gz). See the [release notes](https://bitbucket.org/malb/m4rie/wiki/M4RIE-20111004) for the list of changes.

* **2011/07/14** A new release of M4RI is available for [download](https://bitbucket.org/malb/m4ri/downloads/m4ri-20110715.tar.gz). See the [release notes](https://bitbucket.org/malb/m4ri/wiki/M4RI-20110715) for the list of changes. Also, a new release of M4RIE is also available for [download](https://bitbucket.org/malb/m4ri/downloads/m4rie-20110715.tar.gz). M4RIE now relies on M4RI for cache size and other hardware feature detection.

* **2011/06/10** A new release of M4RI is available for [download](https://bitbucket.org/malb/m4ri/downloads/m4ri-20110613.tar.gz). This version fixes various issues when M4RI is built with OpenMP enabled.

* **2011/06/01** A new release of M4RI is available for [download.](https://bitbucket.org/malb/m4ri/downloads/m4ri-20110601.tar.gz) See the [release notes](https://bitbucket.org/malb/m4ri/wiki/M4RI-20110601) for the list of changes. Also, a new release of M4RIE is also available for [download](https://bitbucket.org/malb/m4ri/downloads/m4rie-20110601.tar.gz). The only changes to M4RIE are to ensure compatibility with M4RI version 20110601 and up.

* **2011/04/13** We now have a [mailinglist](http://groups.google.com/group/m4ri-devel).

* **2010/08/14** A new release of M4RI is available for [download.](https://bitbucket.org/malb/m4ri/downloads/m4ri-20100817.tar.gz) The main changes are improved automatic cache size detection and some clean ups necessary for M4RIE. A first official release of M4RIE is also available for [download](https://bitbucket.org/malb/m4ri/downloads/m4rie-20100817.tar.gz).

* **2010/07/13** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20100701.tar.gz) is available for download. See the [release notes](http://www.bitbucket.org/malb/m4ri/wiki/M4RI-20100701) for details.

* **2009/11/04** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20091101.tar.gz) is available for download. See the [release notes](http://www.bitbucket.org/malb/m4ri/wiki/M4RI-20091101) for details.

* **2009/04/09** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20090409.tar.gz) is available for download. It heavily breaks backward compatibility but supports much bigger matrices than before. See the [release notes](http://www.bitbucket.org/malb/m4ri/wiki/M4RI-20090409) for details.

* **2009/01/05** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20090105.tar.gz) is available for download. It contains new features, performance enhancements and bug fixes. [Release notes](http://www.bitbucket.org/malb/m4ri/wiki/M4RI-20090105) are available in the wiki.

* **2008/11/12** A paper describing our matrix multiplication implementation is available as [pre-print](http://arxiv.org/abs/0811.1714) on the ArXiv. Also, M4RI is being [packaged](https://bugzilla.redhat.com/show_bug.cgi?id=470173) for Fedora Core. Finally, we updated the [peformance](./performance.html) data for GAP and Magma on the Core 2 Duo with improved timings.

* **2008/10/28** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20081029.tar.gz) is available for download. It contains mainly bugfixes but starting with this release triangular solving with matrices (TRSM) is fully supported. Also LUP factorisation (i.e. on full rank matrices) seems to be working now but it is not optimised at all.

* **2008/10/22** The [slides](http://www.informatik.uni-bremen.de/~malb/talks/20081010%20-%20M4RI%20-%20Nancy.pdf) for the [Sage Days 10](http://wiki.sagemath.org/days10) talk about matrix multiplication in the M4RI library are available online.

* **2008/09/22** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20080909.tar.gz) is available. It is identical to the version of M4RI shipped with [Sage](http://www.sagemath.org) 3.1.2 and contains many build fixes for a wide range of platforms. Sage (and thus M4RI) supportes x86 Linux, x86_64 Linux, ia64 Linux, x86 OSX and ppc OSX. M4RI also supports Windows and Solaris 10.

* **2008/08/26** This [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20080826.tar.gz) is a pure bugfix release. Before this bugfix, if the input matrices were very non-square either wrong results or SIGSEGVs could be observed.

* **2008/08/21** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20080821.tar.gz) is available. This release contains Clément Pernet's latest LQUP and TRSM development code. LQUP still lacks a basecase but TRSM should be fairly complete. No attempts were made so far to optimise things. Furthermore, this release contains an improved strategy for choosing $k$ in M4RM which improves performance on the Core2Duo.

* **2008/08/17** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20080817.tar.gz) is available. This release adds a simple memory manager for systems with slow malloc/free syscalls. Also, the initialisation (m4ri_init) and finalisation (m4ri_fini) routines are now called automatically when the library is loaded/unloaded. This is tested with GCC and SunCC but not with MSVC. Matrix elimination got slightly [faster](./performance.html) across plattforms, multi-core support was extended to elimination and improved for multiplication. The README contains instruction how to enable multi-core support. This release does not contain Clément Pernet's latest LQUP patch.

* **2008/06/24** A new [release](https://bitbucket.org/malb/m4ri/downloads/m4ri-20080624.tar.gz) is available. This release uses the libtool [-release mechanism](http://www.gnu.org/software/libtool/manual.html#Versioning) to ensure binary (in)compatibility between releases since - again - the API changed: since the project is quite young do not expect the API to be stable anytime soon. Also the new version [attempts to detect](http://autoconf-archive.cryp.to/ax_cache_size.html) the L1 and L2 cache sizes and uses a Strassen-Winograd cutoff by default such that both source matrices fit in L2 (this is not optimal but a good compromise). This new version has some scratch/experimental code which is the beginning of asymptotically fast LQUP factorisation. Finally, elimination got slightly faster.

* **2008/06/20** Thanks to [Tim Abbott](http://web.mit.edu/tabbott/www/) libM4RI is now in Debian/unstable.

* **2008/06/13** It turns out our [comparison with Magma](./performance.html) on the Core2Duo was strongly biased, since we compared with a version of Magma that was optimised for AMD64 rather than Intel64\. The correct times are given now and we apologise for this mix-up.

* **2008/06/03** This release is a small bugfix release. Matrices are now printed correctly and a bug in mzd_gauss_delayed was reported and fixed by Wael Said and Mohammed Saied.

* **2008/06/01** This release greatly improves the performance of M4RI: the reduction of a given matrix to (reduced) row echelon form. The speed-up over the last release can be as much as ten, we will provide performance data for this in the near future. However, the new implementation still isn't asymptotically fast. Also mzd_transpose is much faster now due to improved data locality.

* **2008/05/21** Today's release fixes a severe bug found by Bill Hart, disables SSE2 on all CPUs except those manufactured by Intel (for performance reasons), improves [performance](./performance.html) on the Core2Duo and introduces a configure switch to enable OpenMP support.

* **2008/05/20** A new release is available with [massive speed improvements](./performance.html) for matrix multiplication. These improvements were discussed and tested in [this thread](http://groups.google.com/group/sage-devel/browse_thread/thread/aa4edc241ca4d6bb) on the [sage-devel](http://groups.google.com/group/sage-devel) mailing list. This release has also experimental and preliminary support for [OpenMP](http://www.openmp.org). To activate it compile with GCC 4.2 and `CFLAGS="-fopenmp -DHAVE_OPENMP“` Note however, that this release is still a developer preview since some automatic tuning is still not implemented, the performance on the Opteron isn't acceptable yet, and the parallel implementation is naive.

* **2008/05/16** Release early, release often. This release fixes the unconditional use of _mm_free even when it is not available.

* **2008/05/15** A new minor release is available which improves performance on Opterons. Also, the website has moved to [http://m4ri.sagemath.org](http://m4ri.sagemath.org/).
