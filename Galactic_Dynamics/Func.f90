module Func
contains

  function f(r, rs1, ms1)
    use Global
    implicit none

    real*8, dimension(3) :: f, r, rs1
    real*8 :: d1, d2, ms1
    d1 = sqrt(dot_product(r-rs1, r-rs1))
    f = ms1/d1**3*(rs1-r)
    return
  end function f

  function cross_product_2d(a, b)
    implicit none

    real*8, dimension(2) :: a, b
    real*8 :: cross_product_2d
    cross_product_2d = a(1)*b(2)-a(2)*b(1)
    return
  end function cross_product_2d

  function cross_product_3d(a, b)
    implicit none

    real*8, dimension(3) :: a, b, cross_product_3d
    cross_product_3d(1) = a(2)*b(3)-a(3)*b(2)
    cross_product_3d(2) = a(3)*b(1)-a(1)*b(3)
    cross_product_3d(3) = a(1)*b(2)-a(2)*b(1)
    return
  end function cross_product_3d

  subroutine RK8(f, r, v, t, h, rs1, ms1)
    implicit none

    real*8, dimension(3) :: k1_1, k1_2, k1_3, k1_4, k1_5, k1_6, k1_7, k1_8, k1_9, k1_10
    real*8, dimension(3) :: k2_1, k2_2, k2_3, k2_4, k2_5, k2_6, k2_7, k2_8, k2_9, k2_10
    real*8, dimension(3) :: dv, dr, r, v, rs1
    real*8 :: t, h, ms1

    interface
      function f(r, rs1, ms1)
        real*8, dimension(3) :: f, r, rs1
        real*8 :: d1, ms1
      end function f
    end interface

    k1_1 = f(r, rs1, ms1)
    k2_1 = v

    k1_2 = f(r + (h*4./27.)*k2_1,&
             rs1, ms1)
    k2_2 = v + (h*4./27.)*k1_1

    k1_3 = f(r + (h/18.)*(k2_1+3*k2_2),&
             rs1, ms1)
    k2_3 = v + (h/18.)*(k1_1+3*k1_2)

    k1_4 = f(r + (h/12.)*(k2_1+3*k2_3),&
             rs1, ms1)
    k2_4 = v + (h/12.)*(k1_1+3*k1_3)

    k1_5 = f(r + (h/8.)*(k2_1+3*k2_4),&
             rs1, ms1)
    k2_5 = v + (h/8.)*(k1_1+3*k1_4)

    k1_6 = f(r + (h/54.)*(13*k2_1-27*k2_3+42*k2_4+8*k2_5),&
             rs1, ms1)
    k2_6 = v + (h/54.)*(13*k1_1-27*k1_3+42*k1_4+8*k1_5)

    k1_7 = f(r + (h/4320.)*(389*k2_1-54*k2_3+966*k2_4-824*k2_5+243*k2_6),&
             rs1, ms1)
    k2_7 = v + (h/4320.)*(389*k1_1-54*k1_3+966*k1_4-824*k1_5+243*k1_6)

    k1_8 = f(r + (h/20.)*(-234*k2_1+81*k2_3-1164*k2_4+656*k2_5-122*k2_6+800*k2_7),&
             rs1, ms1)
    k2_8 = v + (h/20.)*(-234*k1_1+81*k1_3-1164*k1_4+656*k1_5-122*k1_6+800*k1_7)

    k1_9 = f(r + (h/288.)*(-127*k2_1+18*k2_3-678*k2_4+456*k2_5-9*k2_6+576*k2_7+4*k2_8),&
             rs1, ms1)
    k2_9 = v + (h/288.)*(-127*k1_1+18*k1_3-678*k1_4+456*k1_5-9*k1_6+576*k1_7+4*k1_8)

    k1_10 = f(r + (h/820.)*(1481*k2_1-81*k2_3+7104*k2_4-3376*k2_5+72*k2_6-5040*k2_7-60*k2_8+720*k2_9),&
              rs1, ms1)
    k2_10 = v + (h/820.)*(1481*k1_1-81*k1_3+7104*k1_4-3376*k1_5+72*k1_6-5040*k1_7-60*k1_8+720*k1_9)

    dv = h/840.*(41*k1_1+27*k1_4+272*k1_5+27*k1_6+216*k1_7+216*k1_9+41*k1_10)
    dr = h/840.*(41*k2_1+27*k2_4+272*k2_5+27*k2_6+216*k2_7+216*k2_9+41*k2_10)

    v = v + dv
    r = r + dr

  end subroutine RK8

end module Func
