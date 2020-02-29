module Func
contains

  function f(r, rs1, rs2, ms1, ms2)
    use Global
    implicit none

    real*8, dimension(2) :: f, r, rs1, rs2
    real*8 :: d1, d2, ms1, ms2
    d1 = sqrt(dot_product(r-rs1, r-rs1))
    d2 = sqrt(dot_product(r-rs2, r-rs2))
    f = ms1/d1**3*(rs1-r) + ms2/d2**3*(rs2-r)
    return
  end function f

  function cross_product_2d(a, b)
    implicit none

    real*8, dimension(2) :: a, b
    real*8 :: cross_product_2d
    cross_product_2d = a(1)*b(2)-a(2)*b(1)
    return
  end function cross_product_2d

  function radius(mass) !unit - solar radius
    implicit none
    real*8 :: mass, a, b, radius
    a = mass / 1.44
    b = mass / 0.00057
    radius = (0.0114 * (a**(-2. / 3.) - a**(2. / 3.))**(1. / 2.) * &
             (1 + 3.5 * b**(-2. / 3.) + b**(-1.))**(-2. / 3.))
    return
  end function radius

  function thetas(q, x)
    implicit none
    real*8 :: q, x, mu, A, sin, thetas
    mu = q / (1 + q)
    A = mu / abs(-x - 1 + mu)**3 + (1 - mu) / abs(-x + mu)**3
    sin = -sqrt((8. / 9. / A) * (1 - 2 / A + 3 * &
                                  sqrt(1 - 8. / 9. / A)))
    thetas = asin(sin) / 2.
    return
  end function thetas

  function phi(x, q)
    implicit none
    real*8 :: phi, x, q
    phi = -1. / (x - q / (1. + q))**2 + q / (x + 1. / (1. + q))**2 - (1. + q) * x
    return
  end function phi

  function phi_prime(x, q)
    implicit none
    real*8 :: phi_prime, x, q
    phi_prime = -1. - q - (2. * q)/(1. / (1. + q) + x)**3 + 2. / (-q / (1. + q) + x)**3
    return
  end function phi_prime

  subroutine ejection(md, m3, rd, r, vd, v, fd, vth)
    use Global
    implicit none
    real*8 :: md, m3, mdf, fd, radius1
    real*8, dimension(2) :: rd, r, rdf, vd, v, vdf, vth
    common/ejection/ radius1
    !!Mass
    mdf = md - m3
    !!Position
    rdf = (md*rd - m3*r) / mdf
    !!Orbital velocity
    v = v + vth
    vdf = (md*vd - m3*v) / mdf
    !!Spin
    fd = (0.4*md*radius1**2 * fd &
        + md*cross_product_2d(rd, vd) &
        - mdf*cross_product_2d(rdf, vdf) &
        - m3*cross_product_2d(r, v)) &
        / (0.4*mdf*radius1**2)
    md = mdf
    rd = rdf
    vd = vdf

  end subroutine ejection

  subroutine accretion(ma, m3, ra, r, va, v, fa)
    use Global
    implicit none
    real*8 :: ma, m3, maf, fa, radius2
    real*8, dimension(2) :: ra, r, raf, va, v, vaf
    common/accretion/ radius2
    !!Mass
    maf = ma + m3
    !!Position
    raf = (ma*ra + m3*r) / maf
    !!Orbital velocity
    vaf = (ma*va + m3*v) / maf
    !!Spin
    fa = (0.4*ma*radius2**2 * fa &
        + ma*cross_product_2d(ra, va)&
        - maf*cross_product_2d(raf, vaf) &
        + m3*cross_product_2d(r, v)) &
        / (0.4*maf*radius2**2)

    ma = maf
    ra = raf
    va = vaf
    m3 = 0

  end subroutine accretion

  subroutine RK8(f, r, v, t, h, rs1, rs2, ms1, ms2)
    implicit none

    real*8, dimension(2) :: k1_1, k1_2, k1_3, k1_4, k1_5, k1_6, k1_7, k1_8, k1_9, k1_10
    real*8, dimension(2) :: k2_1, k2_2, k2_3, k2_4, k2_5, k2_6, k2_7, k2_8, k2_9, k2_10
    real*8, dimension(2) :: r, v, rs1, rs2
    real*8 :: t, h, ms1, ms2

    interface
      function f(r, rs1, rs2, ms1, ms2)
        real*8, dimension(2) :: f, r, rs1, rs2
        real*8 :: d1, d2, ms1, ms2
      end function f
    end interface

    k1_1 = f(r, rs1, rs2, ms1, ms2)
    k2_1 = v

    k1_2 = f(r + (h*4./27.)*k2_1,&
             rs1, rs2, ms1, ms2)
    k2_2 = v + (h*4./27.)*k1_1

    k1_3 = f(r + (h/18.)*(k2_1+3*k2_2),&
             rs1, rs2, ms1, ms2)
    k2_3 = v + (h/18.)*(k1_1+3*k1_2)

    k1_4 = f(r + (h/12.)*(k2_1+3*k2_3),&
             rs1, rs2, ms1, ms2)
    k2_4 = v + (h/12.)*(k1_1+3*k1_3)

    k1_5 = f(r + (h/8.)*(k2_1+3*k2_4),&
             rs1, rs2, ms1, ms2)
    k2_5 = v + (h/8.)*(k1_1+3*k1_4)

    k1_6 = f(r + (h/54.)*(13*k2_1-27*k2_3+42*k2_4+8*k2_5),&
             rs1, rs2, ms1, ms2)
    k2_6 = v + (h/54.)*(13*k1_1-27*k1_3+42*k1_4+8*k1_5)

    k1_7 = f(r + (h/4320.)*(389*k2_1-54*k2_3+966*k2_4-824*k2_5+243*k2_6),&
             rs1, rs2, ms1, ms2)
    k2_7 = v + (h/4320.)*(389*k1_1-54*k1_3+966*k1_4-824*k1_5+243*k1_6)

    k1_8 = f(r + (h/20.)*(-234*k2_1+81*k2_3-1164*k2_4+656*k2_5-122*k2_6+800*k2_7),&
             rs1, rs2, ms1, ms2)
    k2_8 = v + (h/20.)*(-234*k1_1+81*k1_3-1164*k1_4+656*k1_5-122*k1_6+800*k1_7)

    k1_9 = f(r + (h/288.)*(-127*k2_1+18*k2_3-678*k2_4+456*k2_5-9*k2_6+576*k2_7+4*k2_8),&
             rs1, rs2, ms1, ms2)
    k2_9 = v + (h/288.)*(-127*k1_1+18*k1_3-678*k1_4+456*k1_5-9*k1_6+576*k1_7+4*k1_8)

    k1_10 = f(r + (h/820.)*(1481*k2_1-81*k2_3+7104*k2_4-3376*k2_5+72*k2_6-5040*k2_7-60*k2_8+720*k2_9),&
              rs1, rs2, ms1, ms2)
    k2_10 = v + (h/820.)*(1481*k1_1-81*k1_3+7104*k1_4-3376*k1_5+72*k1_6-5040*k1_7-60*k1_8+720*k1_9)

    v = v + h/840.*(41*k1_1+27*k1_4+272*k1_5+27*k1_6+216*k1_7+216*k1_9+41*k1_10)
    r = r + h/840.*(41*k2_1+27*k2_4+272*k2_5+27*k2_6+216*k2_7+216*k2_9+41*k2_10)

  end subroutine RK8

end module Func
