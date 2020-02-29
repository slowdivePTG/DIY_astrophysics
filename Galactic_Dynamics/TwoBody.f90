program Two_Body

  use Func
  use Global
  implicit none
  !real*8, parameter :: RL = a*(1/(1+q) + x)
  !real*8, parameter :: RL1 = a*(q/(1+q) - x)
  !real*8, parameter :: Rc = (1+q)*RL1**4/a**3
  !Randomly generated
  !Other initial parameters
  real*8 :: theta, thetac
  !model parameters
  real*8 :: t, dt, a
  !Mass
  real*8 :: m1, m2, q
  !velocity & location
  real*8, dimension(2) :: v1, v2, r1, r2, r1_temp, r2_temp
  !angular momentum
  real*8 :: j1, j2, j
  integer:: i

  !output header
    open(unit=71, file="./data/R_V")
    open(unit=72, file="./data/J")

    q = 1e-10

    !Initial condition
    !!Mass & Radius
    m1 = q / (1 + q)
    m2 = 1 / (1 + q)

    !!Position
    r1(1) = -1 / (1 + q)
    r1(2) = 0
    r2(1) = q / (1 + q)
    r2(2) = 0

    !!Orbital velocity
    v1(1) = 0
    v2(2) = 0
    v2(1) = 0
    v2(2) = 0

    j1 = m1*cross_product_2d(r1, v1)
    j2 = m2*cross_product_2d(r2, v2)
    j = j1 + j2

    t = 0d0
    write(71,'(a)') '#Time, R1x, R1y, R2x, R2y, V1x, V1y, V2x, V2y'
    write(71,'(999E22.12)') t, r1, r2, v1, v2
    write(72,'(a)') '#Time, J1, J2, J'
    write(72,'(999E22.12)') t, j1, j2

    i = 0
    thetac = 0 !Stop when one orbit completes
    do
      i = i + 1
      if (i > nstep) exit
      if (thetac > pi*2) exit
      dt = time_interval/(max(sqrt(dot_product(v,v)),1.))
      r1_temp = r1
      r2_temp = r2
      call RK8(f, r1, v1, t, dt, r2_temp, m2)
      call RK8(f, r2, v2, t, dt, r1_temp, m1)
      t = t + dt
      thetac = thetac + asin(cross_product_2d(r1,r1_temp)/&
               sqrt(dot_product(r1,r1)*dot_product(r1_temp,r1_temp)))
      
      j1 = m1*cross_product_2d(r1, v1)
      j2 = m2*cross_product_2d(r2, v2)
      j = j1 + j2
      if (mod(i, plot_interval) == 0) then
        write(71,'(999E22.12)') t, r1, r2, v1, v2
        write(72,'(999E22.12)') t, j1, j2, j
      end if
    end do
      
    close(unit=71)
    close(unit=72)

end program Two_Body
