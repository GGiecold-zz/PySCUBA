# PySCUBA

PySCUBA stands for "Python for Single-cell Clustering Using Bifurcation Analysis." PySCUBA is a novel computational method for extracting lineage relationships from single-cell genomics data, and modeling the dynamic changes associated with cell differentiation. PySCUBA draws techniques from nonlinear dynamics and stochastic differential equation theories, providing a systematic framework for modeling complex processes involving multi-lineage specifications. 

Installation
------------

PySCUBA is written in Python 2.7. It relies on the following libraries and APIs:
* ```matplolib``` (1.4.3 or ulterior version)
* ```numpy``` (1.9.0 or later)
* ```Pillow``` (3.2.0 or later)
* ```PyQt4``` (at least version 4.11.4)
* ```python-igraph``` (version>=0.7.1)
* ```rpy2``` (2.8.1 or later)
* ```scipy``` (0.17.0 or later)
* ```setuptools``` (version>=21.0.0)
* ```sklearn``` (version>=0.17.1)
* ```Wand``` (version 0.4.3 or subsequent)

Most of those dependencies should be automatically resolved during an installation with ```pip```, by
* starting a terminal;
* running ```$ pip install PySCUBA```.

Nonetheless, please ensure those prerequisites have been met before proceeding further. In particular, PyQt might have to be installed manually. In addition, note that you might need to run this using ```sudo``` depending on your Python installation.

Alternatively, you can also get PySCUBA from source by downloading a tarball or zip archive from the corresponding GitHub page. Once you have downloaded and unpacked the source, you can navigate into the root source directory and run:

```$ python setup.py install```

Here be dragons!
----------------

This is still work in progress.

I have been tasked with redesigning into a Python package the Matlab code supporting the analysis presented in the following paper:
Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.

The motivation for this project is twofold:
- to make do without high Matlab licensing fees that were preventing some of our colleagues from using SCUBA and other tools developed at our lab;
- to come up with faster, more efficient way of handling the computations and data processing; a benchmark comparison done on the same machine and on a few datasets found that where SCUBA took about 20 minutes to complete, PySCUBA handles the same task in about 30 seconds.

I have found a few bugs, mathematical inconsistencies and improper or overlooked handling of outlier cases in the Matlab code. Such issues are fully addressed in the present Python package.

Attribution
-----------

If you find PySCUBA useful in your research, please cite its GitHub repository:
Giecold G, PySCUBA, (2016), GitHub repository, https://github.com/GGiecold/PySCUBA

The respective BibTex entry is
```
@misc{Giecold2016,
  author = {Giecold, G.},
  title = {PySCUBA},
  year = {2016},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/GGiecold/PySCUBA}},
  commit = {}
}
```

Besides, please cite the following reference as well:
Marco E, Karp RL, Guo G, Robson P, Hart AH, Trippa L, Yuan GC. (2014) "Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape". Proceedings of the National Acacdemy of Sciences, 111, E5643-E5650.

Licence
-------

Copyright 2016-2021 Gregory Giecold.

PySCUBA is free software made available under the MIT License. For details see the LICENSE file.



