# This is an commentary line in a makefile

# Start of the makefile

object = Global.o RadialInfall.o Func.o

ThreeBody: $(object)

	gfortran -o RadialInfall $(object)

Global.mod: Global.o Global.f90

	gfortran -c Global.f90

Global.o: Global.f90

	gfortran -c Global.f90

Func.mod: Func.o Func.f90

	gfortran -c Global.f90

Func.o: Global.mod Func.f90

	gfortran -c Func.f90

RadialInfall.o: Global.mod Func.mod RadialInfall.f90

	gfortran -c RadialInfall.f90

%.o: %.f90
	fortran -c $<

clean:

	rm *.mod *.o RadialInfall

# End of the makefile
