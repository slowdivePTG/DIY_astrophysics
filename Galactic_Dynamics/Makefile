# This is an commentary line in a makefile

# Start of the makefile

object = Global.o TwoBody.o Func.o Newton.o Init.o

TwoBody: $(object)

	gfortran -o TwoBody $(object)

Global.mod: Global.o Global.f90

	gfortran -c Global.f90

Global.o: Global.f90

	gfortran -c Global.f90

Func.mod: Func.o Func.f90

	gfortran -c Global.f90

Func.o: Global.mod Func.f90

	gfortran -c Func.f90

TwoBody.o: Global.mod Func.mod TwoBody.f90

	gfortran -c TwoBody.f90

%.o: %.f90
	fortran -c $<

clean:

	rm *.mod *.o TwoBody

# End of the makefile
