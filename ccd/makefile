#!/bin/sh
f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -m ueg_ccd -c Precision.f90 Constants.f90 ccd.f90 ccdt.f90
#f2py --fcompiler=gfortran -m ueg_ccd -c Precision.f90 Constants.f90 CCD_dr.f90

#rm *mod
#rm *.o
