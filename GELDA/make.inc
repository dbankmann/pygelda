#
#  GELDA make include file.                                        
#  Version 1.1.1 (SunOS)                                             
#  October 12, 1995
#
#  Tested on a Sun SPARCstation 10 under SunOS 4.1.3
#  with the Fortran SPARCompiler Release 2.0.1 from SunPro. 
#  

#----------------------------------------------------------------------
#  YOU MUST EDIT THE FOLLOWING LINES 
#----------------------------------------------------------------------

#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = gfortran
FCOPTS   = 
NOOPT    = 
LOADER   = gfortran
LOADOPTS = 
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  LaTeX
#
LATEX	     = latex
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
#BLASLIB	= ./blas.a
BLASLIB = -L/usr/site-local/lib -lBlasGnu
#LAPACKLIB = ./lapack.a
LAPACKLIB = -L/usr/site-local/lib -lBlasGnu -lLaPackGnu 
#  
#  The source file for the timing code. Furthermore, the compiler, the options
#  and libraries you need to compile and link with the timing code. 
#
TIMESRC = dsecnd.f
TIMECOMP = gfortran
TIMEOPTS = 
TIMELIB = 





