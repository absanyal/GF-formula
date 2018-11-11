rm -r *.dat
rm -r g_level*
clear
gfortran -c *.f95 #-llapack -lblas
gfortran *.o -llapack -lblas
./a.out