program Three_Body

  use Func
  use Global
  implicit none
  integer :: i
  real*8 :: r0, a
  real*8 :: t, r, v
  real*8 :: dt

  !output header
  open(unit=71, file="R_V")

  r0 = 1.0d0
  a = 1.0d-7

  r = r0
  v = 0d0

  t = 0d0
  write(71,'(a)') '#Time, R, V'
  write(71,'(999E22.12)') t, r, v

  i = 0
  do
    i = i + 1

    !Stop
    if (r<=a) then
      exit
    end if

    dt = time_interval/(max(v, 1.))

    call RK8(f, r, v, t, dt)

    if (mod(i, plot_interval) == 0) then
      write(71,'(999E22.12)') t, r, v
    end if

  end do

  close(unit=71)
  close(unit=72)
end program Three_Body
