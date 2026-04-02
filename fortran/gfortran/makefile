srcdir = ../source

all: exec/MTSphere.exe

source = ${srcdir}/precision.f90 ${srcdir}/SphericalFunctions.f90 ${srcdir}/csphjy.f90 ${srcdir}/FILT_COEF_Q.f90 ${srcdir}/HSSPHERE_MT.f90 \
         ${srcdir}/mcikva.f90 ${srcdir}/MTFields.f90 ${srcdir}/MTSphere.f90\
	 ${srcdir}/SphericalHarmonicAnalysis.f90

exec/MTSphere.exe : $(source)
	gfortran -o exec/MTSphere.exe $(source) liblapack.a libblas.a -fstack-protector
