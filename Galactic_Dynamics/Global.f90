module Global
  implicit none

  real*8, parameter :: Rs = 69550001060.63751 !solor radius (cm)
  real*8, parameter :: Ms = 1.98841586d+33 !solor mass (g)
  real*8, parameter :: G = 6.67384d-08 !gravitational constant (cm**3/(g*s**2))
  real*8, parameter :: pi = 3.141592653589793
  real*8, parameter :: spy = 31557600.0 !seconds per year
  real*8, parameter :: Temp = 20000. !Surface temperature of the donor (K)

  integer, parameter :: nstep = 50000000
  integer, parameter :: plot_interval = 10000
  real*8, parameter :: time_interval = 1d-7

end module Global
