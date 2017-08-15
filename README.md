# JACO

## Introduction

The *Joint Analysis of Cluster Observations* project seeks to provide a distributed C-based code (using combined OpenMP and MPI parallelism) to model the matter distribution in galaxy clusters. The authors are Andy Mahdavi, Seth Siegel, Alison Mansheim, and Jack Sayers. Capabilities include:
* Joint fitting of X-ray, weak lensing, and SZ data
* Markov Chain Monte Carlo determination of confidence intervals on matter parameters
* (upcoming) modeling of triaxial halos.

## Prerequisites

JACO runs on Linux systems. On Ubuntu-based systems, the following command would install the prerequisites (requires system administrator access):

`sudo apt-get install openmpi-bin texinfo zlib1g-dev libcfitsio3-dev libfftw3-dev libgsl0-dev libopenmpi-dev pgplot5 libreadline6-dev`

On other Linux systems, use the native package managers (such as `yum`) to install the prerequisites.

## Installation

To install in the system directory (with administrator access), type
`./configure --with-bolocam`
`make -j 4`
`sudo make install`

If you do not have administrative privileges, you can install on your home folder:
`./configure --with-bolocam --prefix=$HOME`
`make -j 4`
`make install`

## Incorporating Seth Siegel's changes

The following files were changed by Seth Siegel, but the differences are too complicated to be managed by GitHub's automatic merge system. Whoever incorporates the changes in these files into the main code should grab relevant portions bookended by ``Siegel Begin`` and ``Siegel End`` (or ``//Siegel``) and manually merge them into the main files. Note that only these changes, and not other changes, should be merged, since outside those specifically marked changes, the codebase has evolved as well.

```
hrothgar/sciutils/readxray_SETH.c
hrothgar/hrothgar_SETH.c
jaco_SETH.c
jaco_SETH.h
sherpa_SETH.c
standalone_SETH.c
standalone_SETH.h
```

## Velocity dispersion modeling

The file `veldisp.c` contains Alison Mansheim's code for modeling velocity dispersions jointly with other code. This code has fallen out of sync with the rest of the code base and while it runs, it will segfault (crash). Some minor updating should be required to get it to work.
