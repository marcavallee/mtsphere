MTSphere Fortran readme.

MTSphere executable has been compiled and executed on 12th Gen Intel Core i7-12700H, 2.3 GHz, with 16 Go RAM memory and a 64 bits system and Windows 11. The computation of the Sounding model (73 frequenciecs) took 44 seconds with the Intel version and 23 seconds if compiled with gfortran. Source code is the same for the compilations.

The input file is named MTSphereMT.cfl. It is an ASCII file with the following lines.

1. Comment line.

2. NLYR, NTERMS, NX, NY, where
	a) NLYR : Layer number
	b) NTERMS : Number of spherical harmonics computed ( 6 recommended )
	c) NX : number of stations along the X direction
	d) NY : number of stations along the Y direction
	If NX and NY are both positive, station locations are interpolated from XMIN and XMAX and YMIN and YMAX. Otherwise, if either NX or NY is negative, the user
	must provide the station locations.

3. NF, MINFREQ, MAXFREQ, LOGARITHMIC, where
	a) NF : Frequency number
	b) MINFREQ : Minimum frequency
	c) MAXFREQ : Maximum frequency
	d) LOGARITHMIC: 0 (linear) or 1 (logarithmic) frequency spacing

4. RADIUS, DEPTH, SPHRES, where
	a) RADIUS : radius of the sphere
	b) DEPTH : Depth of the center of the sphere

5.x NLYR-1 lines where
	a) RES(I), THK(I) : resistivity and thickness of each layer

6. one line with RES(NLYR) : bottom half-space resistivity

If NX or NY is negative:

7. X(1) ... X(NX) : X coordinates of the MT stations

8. Y(1) ... Y(NY) : Y coordinates of the MT stations

Otherwise:

7. XMIN XMAX

8. YMIN YMAX

Results are output as ASCII files.

MTSphere was compiled with Classic Intel Compiler 2021.7.1. imbedded in Visual Studio  and Cygwin gfortran for Windows. A Visual Studio Project is provided. The following instructions are provided for compilation with gfortran.

1) gfortran must be available.

2) blas libraries must be installed. They are provided as is with the distribution.

3) With g fortran, the compilation is simply done by typing Make -f makefile.

July 2025.


