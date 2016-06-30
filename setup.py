#!/usr/bin/env python


# PySCUBA/setup.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com, ggiecold@jimmy.harvard.edu


from codecs import open
from os import path

from setuptools import setup


exec(open(path.join('src', 'PySCUBA', '__version__.py')).read())


here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README'), encoding = 'utf-8') as f:
    long_description = f.read()
    

setup(name = 'PySCUBA',
      version = __version__,
      
      description = "Python for Single-cell Clustering Using Bifurcation Analysis",
      long_description = long_description,
                    
      url = 'https://github.com/GGiecold/PySCUBA',
      download_url = 'https://github.com/GGiecold/PySCUBA',
      
      author = 'Gregory Giecold',
      author_email = 'g.giecold@gmail.com',
      maintainer = 'Gregory Giecold',
      maintainer_email = 'ggiecold@jimmy.harvard.edu',
      
      license = 'MIT License',
      
      platforms = ('Any',),
      install_requires = ['matplotlib>=1.4.3', 'numpy>=1.9.0', 'Pillow>=3.2.0', 
                          'python-igraph', 'rpy2>=2.8.1', 'scipy>=0.17.0', 
                          'setuptools', 'sklearn', 'Wand>=0.4.3'],
                          
      classifiers = ['Development Status :: 4 - Beta',
                     'Environment :: Console',
                     'Intended Audience :: End Users/Desktop',
                     'Intended Audience :: Developers',
                     'Intended Audience :: Healthcare Industry',
                     'Intended Audience :: Science/Research',          
                     'License :: OSI Approved :: MIT License',
                     'Natural Language :: English',
                     'Operating System :: MacOS :: MacOS X',
                     'Operating System :: Microsoft',
                     'Operating System :: POSIX',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     'Topic :: Scientific/Engineering :: Visualization',
                     'Topic :: Scientific/Engineering :: Mathematics',
                     'Topic :: Software Development :: User Interfaces', ],
                   
      packages = ['PySCUBA'],
      package_dir = {'PySCUBA': 'src/PySCUBA'},

      keywords = "bioinformatics biology clustering cytometry gap-statistics "
                 "genomics machine-learning pattern-recognition PCR principal-curve "
                 "qPCR RNASeq single-cell time-series unsupervised-learning",
                 
      entry_points = {
          'console_scripts': ['PySCUBA = PySCUBA.__main__:main'],
          }    
)


