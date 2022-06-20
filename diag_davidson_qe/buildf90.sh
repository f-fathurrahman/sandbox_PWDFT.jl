gfortran -cpp -fPIC -shared -c my_cegterg_v61.f90
gfortran -cpp -fPIC -shared -c cdiaghg.f90
gfortran -fPIC -shared my_cegterg_v61.o cdiaghg.o -o libdavidson.so -lblas -llapack

