# Further Reading

## Matrix Multiplication

The matrix multiplication machinery in the M4RI library is described in [ABH08].

**[ABH08]** Martin Albrecht, Gregory Bard, William Hart. _Algorithm 898: Efficient Multiplication of Dense Matrices over GF(2)_. ACM Transactions on Mathematical Software (TOMS), 2010\. pre-print available at [http://arxiv.org/abs/0811.1714](http://arxiv.org/abs/0811.1714)

### “Method of the Four Russians” Multiplication

This algorithm was created by [ADKF70] for the transitive closure of a graph, then [AHU72] extended it to matrix multiplication over the boolean semiring. A discussion leading to many improvements to the implementation of that algorithm can be followed in [AH08]. A technical report discussing these implementation issues is [ABH08].

**[ADKF70]** V. Arlazarov, E. Dinic, M. Kronrod, and I. Faradzev. _On economical construction of the transitive closure of a directed graph_. Dokl. Akad. Nauk., 194(11), 1970\. (in Russian), English Translation in Soviet Math Dokl.

**[AHU72]** A. Aho, J. Hopcroft, and J. Ullman. _The Design and Analysis of Computer Algorithms_. Addison-Wesley, 1974.

**[AH08]** Martin Albrecht, Bill Hart et al. _slightly OT: new M4RI library_. [sage-devel] mailing list, May 2008\. [http://groups.google.com/group/sage-devel/...](http://groups.google.com/group/sage-devel/browse_thread/thread/aa4edc241ca4d6bb)

### Strassen Multiplication

The general formula is published in [Str69], the strategy for dealing with extra rows and columns was taken from [BHS08] and the operation schedule from [B08] and [DP08]. An overview of the implementation in the M4RI library is [ABH08].

**[Str69]** Volker Strassen. _Gaussian elimination is not optimal_. Numerische Mathematik, 13:354–356, 1969.

**[B08]** Marco Bodrato. _A Strassen-like matrix multiplication suited for squaring and higher power computation_. Proceedings of the 2010 International Symposium on Symbolic and Algebraic Computation, 2010\. [http://marco.bodrato.it/papers/Bodrato2010-StrassenLikeMatrixMultiplicationForSquares.pdf](http://marco.bodrato.it/papers/Bodrato2010-StrassenLikeMatrixMultiplicationForSquares.pdf)

**[BHS08]** Robert Bradshaw, David Harvey and William Stein. _strassen_window_multiply_c_. strassen.pyx, Sage 3.0, 2008\. [http://www.sagemath.org](http://www.sagemath.org)

**[DP08]** Jean-Guillaume Dumas and Clèment Pernet. _Memory efficient scheduling of Strassen-Winograd's matrix multiplication algorithm_. [arXiv:0707.2347](http://arxiv.org/abs/0707.2347), 2008.

## Echelon Forms

M4RI provides two algorithms for computing echelon forms: M4RI and PLE decomposition.

### PLE Decomposition

The PLE decomposition and its relationship to other decompositions such as PLUQ is described in [JPS11]. A draft describing our asymptotically fast PLE decomposition implementation and our implementation of block-iterative PLE decomposition is available as [ABP11]. A pre-print of a previous version of this document discussing the Method of Many People Factorisation (MMPF) is available as [AP10].

**[ABP11]** Martin Albrecht, Gregory Bard and Clément Pernet. _Efficient Dense Gaussian Elimination over the Finite Field with Two Elements_. [arXiv:1111.6549v1](http://arxiv.org/abs/1111.6549), 2011.

**[AP10]** Martin Albrecht and Clément Pernet. _Efficient Decomposition of Dense Matrices over GF(2)_. [arxiv.org:1006.1744](http://arxiv.org/abs/1006.1744), 2011.

**[JPS11]** Claude-Pierre Jeannerod, Clément Pernet and Arne Storjohann. _Rank-profile revealing Gaussian elimination and the CUP matrix decomposition_. [arxiv.org:1112.571](http://arxiv.org/abs/1112.5717), 2011.

### “Method of the Four Russians” Inversion

The algorithm was first published in [Bar06] and is also covered more detailed in [Bar07]. Some implementation details can be found in [AP10].

**[Bar06]** Gregory V. Bard. _Accelerating Cryptanalysis with the Method of Four Russians_. Cryptology ePrint Archive, [Report 2006/251](http://eprint.iacr.org/2006/251), 2006.

**[Bar07]** Gregory V. Bard. _Algorithms for Solving Linear and Polynomial Systems of Equations over Finite Fields with Applications to Cryptanalysis_. Phd thesis, University of Maryland, 2007.

# Other Publications Referencing M4RI

**[S13]** Willem Sonke. _Reconstructing a level-1-network from quartets_. [Bachelor project](http://alexandria.tue.nl/extra1/afstversl/wsk-i/sonke2013.pdf) at Eindhoven University of Technology, 2013.

**[HS13]** Yann Hamdaoui and Nicolas Sendrier. _A Non Asymptotic Analysis of Information Set Decoding_. [IACR ePrint Report 2013/162](https://eprint.iacr.org/2013/162.pdf), 2013.

**[BDEZ12]** Razvan Barbulescu, Jeremie Detrey, Nicolas Estibals and Paul Zimmermann. _Finding Optimal Formulae for Bilinear Maps_. Proceedings of Arithmetic of Finite Fields. LNCS 7369/2102, pages 168-186, Springer 2012.

**[BR12]** Enrico Bertolazzi and Anna Rimoldi. _Fast matrix decomposition in F2_. [arxiv.org:1209.5198](http://arxiv.org/abs/1209.5198), 2012.

**[DP12]** Jean-Guillaume Dumas and Clément Pernet. _Computational linear algebra over finite fields_. Technical report, 2012.

**[TL12]** David B. Thomas and Wayne Luk. _The LUT-SR Family of Uniform Random Number Generators for FPGA Architectures_. IEEE Transactions on Very Large Scale Integration Systems, 2012.

**[Ull12]** Ehsan Ullah. _New Techniques for Polynomial System Solving_. Phd thesis at Universität Passau, Germany, 2102.

**[PTBW11]** Albrecht Petzoldt, Enrico Thomae, Stanislav Bulygin and Christopher Wolf. _Small Public Keys and Fast Verification for Multivariate Quadratic Public Key Systems_. Proceedings of CHES 2011\. LNCS Volume 6917/2011, pages 475-490, Springer 2011.

**[Moh11a]** Mohamed Saied Emam Mohamed. _Improved Strategies for Solving Multivariate Polynomial Equation Systems over Finite Fields_. PhD thesis at TU Darmstadt, 2011.

**[Moh11b]** Wael Said Abdelmageed Mohamed. _Improvements for the XL Algorithm with Applications to Algebraic Cryptanalysis_. PhD thesis at TU Darmstadt. 2011.

**[PFJWA11]** Agoston Petz, Chien-Liang Fok, Christine Julien, Brenton Walker and Calvin Ardi. _Network Coded Routing in Delay Tolerant Networks: An Experience Report_. Proceedings of the 3rd Extreme Conference on Communication (ExtremeCom), 2011.

**[WAPRJ11]** Brenton Walker, Calvin Ardi, Agoston Petz, Jung Ryu and Christine Julien. _Experiments on the spatial distribution of network code diversity in segmented DTNs_. CHANTS '11 Proceedings of the 6th ACM workshop on Challenged networks, 2011.

**[Alb10]** Martin Albrecht. _Algorithmic Algebraic Techniques and their Application to Block Cipher Cryptanalysis_. Phd thesis at University of London, 2010

**[Bri10]** Michael Brickenstein. _Boolean Gröbner bases - Theory, Algorithms and Applications_. Phd thesis at TU Kaiserslautern, Germany, 2010

**[BBDMW10]** Johannes Buchmann, Stanislav Bulygin, Jintai Ding, Wael Said Abd Elmageed Mohamed and Fabian Werner. _Practical Algebraic Cryptanalysis for Dragon-Based Cryptosystems_. Proceedings of 9th International Conference, CANS 2010\. 2010.

**[BCDM10]** Johannes Buchmann, Daniel Cabarcas, Jintai Ding and Mohamed Saied Emam Mohamed. _Flexible Partial Enlargement to Accelerate Gröbner Basis Computation over F2_. Progress in Cryptology – AFRICACRYPT 2010, 2010\. [http://www.springerlink.com/content/7736225qv8741j3v/](http://www.springerlink.com/content/7736225qv8741j3v/)

**[Dem10]** Denise Demirel. _Effizientes Lösen linearer Gleichungssysteme über GF(2) mit GPUs_. Diplomarbeit at TU Darmstadt, 2010.

**[FL10]** Jean-Charles Faugère and Sylvain Lachartre. _Parallel Gaussian elimination for Gröbner bases computations in finite fields_. Proceedings of the 4th International Workshop on Parallel and Symbolic Computation, 2010\. [http://portal.acm.org/citation.cfm?id=1837210.1837225](http://portal.acm.org/citation.cfm?id=1837210.1837225)

**[ACPS09]** Benny Applebaum, David Cash, Chris Peikert and Amit Sahai. _Fast Cryptographic Primitives and Circular-Secure Encryption Based on Hard Learning Problems_. Advances in Cryptology - CRYPTO 2009, 2009\. [http://www.springerlink.com/content/mg304504208wq584/](http://www.springerlink.com/content/mg304504208wq584/)

**[Bar09]** Gregory V. Bard. _Algebraic Cryptanalysis_. Springer Verlag, 2009.

**[BW09]** Nikhil Bansal and Ryan Williams. _Regularity Lemmas and Combinatorial Algorithms_. 50th Annual IEEE Symposium on Foundations of Computer Science, 2009\. [http://www.computer.org/portal/web/csdl/doi/10.1109/FOCS.2009.76](http://www.computer.org/portal/web/csdl/doi/10.1109/FOCS.2009.76)

**[MCDBB09]** Mohamed Saied Emam Mohamed, Daniel Cabarcas, Jintai Ding, Johannes Buchmann and Stanislav Bulygin. _MXL3: An Efficient Algorithm for Computing Gröbner Bases of Zero-Dimensional Ideals_. Information, Security and Cryptology – ICISC 2009, 2009\. [http://www.springerlink.com/content/v0p5729851q404j7/](http://www.springerlink.com/content/v0p5729851q404j7/)

**[BB08]** Stanislav Bulygin and Michael Brickenstein. _Obtaining and Solving Systems of Equations in Key Variables only for the Small Variants of AES_. Cryptology ePrint Archive, Report 2008/435, 2008\. [http://eprint.iacr.org/2008/435](http://eprint.iacr.org/2008/435)

**[MDB08]** Mohamed Saied Emam Mohamed, Jintai Ding and Johannes Buchmann. _Algebraic Cryptanalysis of MQQ Public Key Cryptosystem by MutantXL_. Cryptology ePrint Archive, Report 2008/451, 2008\. [http://eprint.iacr.org/2008/451](http://eprint.iacr.org/2008/451)

**[MMDB08]** Mohamed Saied Emam Mohamed, Wael Said Abd Elmageed Mohamed, Jintai Ding and Johannes Buchmann. _MXL2: Solving Polynomial Equations over GF(2) Using an Improved Mutant Strategy_. Proceedings of Post-Quantum Cryptography 2008, 2008 [http://www.cdc.informatik.tu-darmstadt.de/reports/reports/MXL2.pdf](http://www.cdc.informatik.tu-darmstadt.de/reports/reports/MXL2.pdf)

**[TL09]** Jérémie Tharaud and Raphaël Laurent. _Linear algebra over the field with two elements using GPUs_. [Technical report](http://moais.imag.fr/membres/jean-louis.roch/perso_html/transfert/2009-06-19-IntensiveProjects-M1-SCCI-Reports/tharaud_laurent_article.pdf), 2009.

Google Scholar [search](https://scholar.google.com/scholar?q=%28%22M4RI%22+AND+%22library%22%29) for references to the M4RI librariy.

# How To Cite Us

If you use our libraries in a non-trivial part of your research please consider citing them as follows:

<pre>@manual{M4RI,
    key          = "M4RI",
    author       = "Martin Albrecht and Gregory Bard",
    organization = "The M4RI~Team",
    title        = "{The M4RI Library -- Version 20121224}",
    year         = 2012,
    url          = "\url{http://m4ri.sagemath.org}",
}
</pre>
and cite the appropriate publications mentioned above.

