# PySCUBA
Python for Single-cell Clustering Using Bifurcation Analysis.

Here be dragons!
----------------

This is still work in progress.

I have been tasked with redesigning into a Python package the Matlab code supporting the analysis presented in the following paper:
Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.

The motivation for this project is twofold:
- to make do without high Matlab licensing fees that were preventing some of our colleagues from using SCUBA and other tools developed at our lab;
- to come up with faster, more efficient way of handling the computations and data processing; a benchmark comparison done on the same machine and on a few datasets found that where SCUBA took about 20 minutes to complete, PySCUBA handles the same task in about 30 seconds.

I have found a few bugs, mathematical inconsistencies and improper or overlooked handling of outlier cases in the Matlab code. Such issues are fully addressed in the present Python package.

Reference
---------

Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.



