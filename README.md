# mtsphere

Programs in Fortran and Python have been developed to simulate the magnetotelluric response of a sphere in a layered earth.  Originally designed for airborne  electromagnetic (EM) exploration, where the source is a dipole (see Vallée and Moussaoui, 2023, Exploration Geophysics 54(4), 362-375), it has been adapted to plane waves and coded in Fortran and Python. Based on developments in the evaluation of derivatives with expansions of spherical and cylindrical functions, this program is open source and can be used to validate numerical models or to better understand the MT response of buried conductors in a layered environment. This approach has been submitted to Geophysics for publication.



Source codes and executables are provided for Fortran and Python. Windows Fortran executables are provided in Fortran/\*/exec directories for Intel and gfortran compiler. Python executable was created using pyInstaller. Of course, the best performance is achieved with Fortran.

