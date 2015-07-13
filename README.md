# pca-utils

A set of utilities to quantify separations between classes of observations
in the scores of principal component analysis (PCA), partial least squares
(PLS) and orthogonal projections to latent structures (OPLS) models. The
original peer-reviewed publication introducing **pca-utils** in published
in:

> Worley, B., Powers, R., _Utilities for Quantifying Separation in PCA/PLS-DA
> Scores_, Analytical Biochemistry, 2013, 433(2): 102-104.

## Introduction

Modeling algorithms like PCA, PLS and OPLS project a set of _K_-variate
observations into a low-dimensional "latent" space. In this space, the
original observations are represented by points (scores), and distances
between the points are related to the original distances between the
high-dimensional observations.

Some obvious questions that arise when discriminating between classes are:

1. How far apart are classes in scores space?
2. Are the scores-space separations significant?
3. Is there a higher pattern to the separations?

The **pca-utils** project provides a set of executables that answer these
questions. The executables allow for generating dendrograms, distance
matrices, class ellipses and ellipsoids based on a set of scores.

## Installing

This software is highly portable. As long as you have a recent enough glibc,
pca-utils should compile and install without incident. You can build and
install **pca-utils** as follows:

> git clone git://github.com/geekysuavo/pca-utils.git

> cd pca-utils

> make

> sudo make install

By default, the **pca-utils** binaries will be installed into _/usr/bin_,
and their manual pages will be installed into _/usr/share/man/man1_. If you
need to install somewhere else, you'll need to modify the Makefile.

For more information on how to use the **pca-utils** binaries, please consult
the manual pages.

## Licensing

The **pca-utils** project is released under the [GNU GPL 3.0] (
http://www.gnu.org/licenses/gpl-3.0.html)

The idea is to advance the state of the art in the field, so as long as you
adhere to the requirements of the above license, just have fun with the code!

*~ Brad.*

