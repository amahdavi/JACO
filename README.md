# JACO

## Introduction

The *Joint Analysis of Cluster Observations* project seeks to provide a distributed C-based code (using combined OpenMP and MPI parallelism) to model the matter distribution in galaxy clusters. The authors are Andy Mahdavi, Seth Siegel, Alison Mansheim, and Jack Sayers. Capabilities include:
* Joint fitting of X-ray, weak lensing, and SZ data
* Markov Chain Monte Carlo determination of confidence intervals on matter parameters
* (upcoming) modeling of triaxial halos.

## Prerequisites

JACO runs on Linux systems. On an Ubuntu-based systems, the following command would install the prerequisites (requires system administrator access):

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
