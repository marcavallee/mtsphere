Description

Introduction

pyMTSphere is a Python program for computing the magnetotelluric response of a sphere in a layered earth. It follows
the structure of MTSphere Fortran program. On the same computer, the performance is better in Fortran.

Algorithm

The algorithm followed is a variant of the algorithm presented by Vallée and Moussaoui (2023). The primary plane-wave
field is developed as a sum of spherical functions following the development of Ward and Hohmann (1987), p. 291-296.
Singular spherical potentials are transformed to cylindrical potentials and layered earth correction are computed
following the developments presented by Vallée and Moussaoui (2023, 2025). Spherical harmonic analysis is required for
computing the layered earth correction. Secondary fields are computed using Hankel transform provided by the python module empymod (Werthmuller 2017).

Programs

The Python functions required for computing the electromagnetic response are provided with this file and include:
. mtplot.py
. layerearthfunctions.py
. matsphere3d.py
. halfspacesphere.py
. hsphere_mt.py
. layerearthspherecorrection.py
. sphericalfunctions.py.

The structure of the various functions follows the Fortran version structure.

References

Vallée, Marc A. and Mouhamed Moussaoui. 2023. Modelling the electromagnetic response of a sphere located in a layered earth.
Exploration Geophysics, 54(4), 362–375.

Vallée, Marc A. and Mouhamed Moussaoui. 2025. Program simulating semi-analytically the magnetotelluric response of a sphere in a layered earth. Submitted to Geophysics.

Ward, S.H., and Hohmann, G.W., 1987, Electromagnetic theory for geophysical applications.
in Electromagnetic methods in applied geophysics, Vol. 1, Theory, ed. M.N. Nabighian. SEG Books.

Werthmuller, D. 2017. An open-source full 3D electromagnetic modeler for 1D VTI media in Python: empymod. Geophysics 82:6 WB9-WB19.