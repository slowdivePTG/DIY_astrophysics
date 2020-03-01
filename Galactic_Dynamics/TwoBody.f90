program Two_Body

  use Func
  use Global
  implicit none
  !initial velocity
  real*8 :: v10, v20, theta1, phi1, theta2, phi2
  real*8 :: thetac
  !model parameters
  real*8 :: t, dt, a
  !Mass
  real*8 :: m1, m2, q
  !velocity & location
  real*8, dimension(3) :: v1, v2, r1, r2, r1_temp, r2_temp, j1, j2
  !angular momentum
  integer:: i

  !output header
  open(unit=71, file="./data_task2/R_V")
  open(unit=72, file="./data_task2/J")
  open(unit=73, file="./data_task2/Para")

  q = 2.d-1
  v10 = 1.d0
  theta1 = pi/3.
  phi1 = pi/3.
  v20 = 1.d-1
  theta2 = pi/4.
  phi2 = -pi/12.

  !Initial condition
  !!Mass & Radius
  m1 = q / (1. + q)
  m2 = 1. / (1. + q)

  !!Position
  r1(1) = -1. / (1. + q)
  r1(2) = 0.
  r1(3) = 0.
  r2(1) = q / (1. + q)
  r2(2) = 0.
  r2(3) = 0.

  !!Orbital velocity
  v1(1) = v10*sin(theta1)*cos(phi1)
  v1(2) = v10*sin(theta1)*sin(phi1)
  v1(3) = v10*cos(theta1)
  v2(1) = v20*sin(theta2)*cos(phi2)
  v2(2) = v20*sin(theta2)*sin(phi2)
  v2(3) = v20*cos(theta2)

  j1 = m1*cross_product_3d(r1, v1)
  j2 = m2*cross_product_3d(r2, v2)

  t = 0d0
  write(71,'(a)') '#Time, R1x, R1y, R2x, R2y, V1x, V1y, V2x, V2y'
  write(71,'(999E22.12)') t, r1, r2, v1, v2
  write(72,'(a)') '#Time, J1, J2, J'
  write(72,'(999E22.12)') t, j1, j2
  write(73,'(a)') '#q, v10, v20, theta1, phi1, theta2, phi2'
  write(73,'(999E22.12)') q, v10, v20, theta1, phi1, theta2, phi2

  i = 0
  thetac = 0 !Stop when one orbit completes
  do
    i = i + 1
    if (i > nstep) then
      print *, 'a'
      exit
    end if
    if (thetac > pi*2) then
      print *, 'b'
      exit
    end if
    if (sqrt(dot_product(r1-r2, r1-r2)) < 1.e-5) then
      print *, 'c', r1, r2
      exit
    end if
    dt = time_interval &
      / (max(max(sqrt(dot_product(v1,v1)),sqrt(dot_product(v2,v2))), 1.))
    r1_temp = r1
    r2_temp = r2
    call RK8(f, r1, v1, t, dt, r2_temp, m2)
    call RK8(f, r2, v2, t, dt, r1_temp, m1)
    t = t + dt
    thetac = thetac + asin(cross_product_2d(r1,r1_temp)/&
             sqrt(dot_product(r1,r1)*dot_product(r1_temp,r1_temp)))
    
    j1 = m1*cross_product_3d(r1, v1)
    j2 = m2*cross_product_3d(r2, v2)
    if (mod(i, plot_interval) == 0) then
      write(71,'(999E22.12)') t, r1, r2, v1, v2
      write(72,'(999E22.12)') t, j1, j2
    end if
  end do
    
  close(unit=71)
  close(unit=72)
  close(unit=73)

end program Two_Body
