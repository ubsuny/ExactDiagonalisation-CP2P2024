objects = main.o
#parameters.o matrixops.o analyze.o main.o
flags = 
#-fcheck=all #this flag doesn't seem to work with mpifort!
compiler = mpif90
#mpifort#gfortran#
LAPACK = -L/network/rit/misc/software/lapack-3.8.0/
LBLAS = -L/network/rit/misc/software/lapack-3.8.0/BLAS/

main: $(objects)
	$(compiler) $(flags) $(objects) -o main.exe $(LBLAS) -lblas $(LAPACK) -llapack
# $(LBLAS) -lblas

main.o: main.f90
	$(compiler) -c $(flags) main.f90

clean:
	rm *.exe *.o *.mod

extra_clean:
	rm *.dat
