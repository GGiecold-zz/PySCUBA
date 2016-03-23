# PySCUBA
Python for Single-cell Clustering Using Bifurcation Analysis

Here be dragons!
----------------

This is still work in progress. 
I have been tasked with redesigning into a Python package the Matlab code supporting the analysis presented in the following paper:
Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.
The motivation for this project is manifold:
- to make away with high Matlab license fees;
- to come up with faster, more efficient way of handling the computations and data processing; a benchmark comparison done on the same machine and on a few datasets found that where SCUBA took about 20 minutes to complete, PySCUBA handles the same task in about 30 seconds.
I have found in the process a few bugs and mathematical inconsistencies that have been fully addressed in the present PySCUBA pacakge.

Reference
---------

Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.



