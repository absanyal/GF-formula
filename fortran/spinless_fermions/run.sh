rm -r *.dat
clear
gfortran -c *.f95 #-llapack -lblas
gfortran *.o -llapack -lblas
./a.out