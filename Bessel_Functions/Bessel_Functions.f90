!-----------------------------------------------------------------------!
!                                                                       !
!                          Bessel_Functions                             !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Bessel_Functions

implicit none

REAL            :: pi = 3.1415926535897932d0
REAL            :: seed
REAL, PARAMETER :: smallest_number = tiny(seed)
REAL, PARAMETER :: biggest_number = huge(seed)

interface BesselJ
  module procedure BesselJ_IR, BesselJ_RR, BesselJ_IC, BesselJ_RC, BesselJ_CR, BesselJ_CC
end interface BesselJ

interface BesselY
  module procedure BesselY_IR, BesselY_RR, BesselY_IC, BesselY_RC, BesselY_CR, BesselY_CC
end interface BesselY

interface BesselI
  module procedure BesselI_IR, BesselI_RR, BesselI_IC, BesselI_RC, BesselI_CR, BesselI_CC
end interface BesselI

interface BesselK
  module procedure BesselK_IR, BesselK_RR, BesselK_IC, BesselK_RC, BesselK_CR, BesselK_CC
end interface BesselK

contains

!-----------------------------------------------------------------------------!
!                                   BesselJ                                   !
!                                                                             !
! This set of functions returns the value of the Bessel function of the first !
! kind J_n(x) for real as well as integer orders n and for real as well as    !
! complex arguments x.                                                        !
!                                                                             !
! Inputs:                                                                     !
!    n: Order of the Bessel function.                                         !
!    x: Argument of the Bessel function.                                      !
!                                                                             !
! Outputs:                                                                    !
!    BesselJ: Value of the Bessel function J_n(x).                            !
!                                                                             !
!-----------------------------------------------------------------------------!

elemental function BesselJ_IR(n,x)
integer, intent(in) :: n
real,    intent(in) :: x
real                :: BesselJ_IR

integer :: nm
real    :: BJ, DJ, BY, DY

call JYNB(n, x, nm, BJ, DJ, BY, DY)

BesselJ_IR = BJ

end function BesselJ_IR

!-----------------------------------------------------------------------------!

elemental function BesselJ_RR(v,x)
real,    intent(in) :: v
real,    intent(in) :: x
real                :: BesselJ_RR

real    :: vm, BJ, DJ, BY, DY

call JYV(v, x, vm, BJ, DJ, BY, DY)

BesselJ_RR = BJ

end function BesselJ_RR

!-----------------------------------------------------------------------------!

elemental function BesselJ_IC(n,z)
integer, intent(in) :: n
complex, intent(in) :: z
complex             :: BesselJ_IC

integer :: nm
complex :: CBJ, CDJ, CBY, CDY

call CJYNB(n, z, nm, CBJ, CDJ, CBY, CDY)

BesselJ_IC = CBJ

end function BesselJ_IC

!-----------------------------------------------------------------------------!

elemental function BesselJ_RC(v,z)
real,    intent(in) :: v
complex, intent(in) :: z
complex             :: BesselJ_RC

real    :: vm
complex :: CBJ, CDJ, CBY, CDY

call CJYVB(v, z, vm, CBJ, CDJ, CBY, CDY)

BesselJ_RC = CBJ

end function BesselJ_RC

!-----------------------------------------------------------------------------!

elemental function BesselJ_CR(zz,v)
complex, intent(in) :: zz
real,    intent(in) :: v
complex             :: BesselJ_CR

complex :: sol, z
z=CMPLX(v,0d0)

call cbsslj(zz,z,sol)

BesselJ_CR = sol

end function BesselJ_CR

!-----------------------------------------------------------------------------!

elemental function BesselJ_CC(zz,z)
complex, intent(in) :: zz
complex, intent(in) :: z  
complex             :: BesselJ_CC

complex :: sol

call cbsslj(zz,z, sol)

BesselJ_CC = sol

end function BesselJ_CC


!-----------------------------------------------------------------------------!
!                                   BesselY                                   !
!                                                                             !
! This set of functions returns the value of the Bessel function of the se-   !
! cond kind Y_n(x) for real as well as integer orders n and for real as well  !
! as complex arguments x.                                                     !
!                                                                             !
! Inputs:                                                                     !
!    n: Order of the Bessel function.                                         !
!    x: Argument of the Bessel function.                                      !
!                                                                             !
! Outputs:                                                                    !
!    BesselY: Value of the Bessel function Y_n(x).                            !
!                                                                             !
!-----------------------------------------------------------------------------!

elemental function BesselY_IR(n,x)
integer, intent(in) :: n
real,    intent(in) :: x
real                :: BesselY_IR

integer :: nm
real    :: BJ, DJ, BY, DY

call JYNB(n, x, nm, BJ, DJ, BY, DY)

BesselY_IR = BY

end function BesselY_IR

!-----------------------------------------------------------------------------!

elemental function BesselY_RR(v,x)
real,    intent(in) :: v
real,    intent(in) :: x
real                :: BesselY_RR

real    :: vm, BJ, DJ, BY, DY

call JYV(v, x, vm, BJ, DJ, BY, DY)

BesselY_RR = BY

end function BesselY_RR

!-----------------------------------------------------------------------------!

elemental function BesselY_IC(n,z)
integer, intent(in) :: n
complex, intent(in) :: z
complex             :: BesselY_IC

integer :: nm
complex :: CBJ, CDJ, CBY, CDY

call CJYNB(n, z, nm, CBJ, CDJ, CBY, CDY)

BesselY_IC = CBY

end function BesselY_IC

!-----------------------------------------------------------------------------!

elemental function BesselY_RC(v,z)
real,    intent(in) :: v
complex, intent(in) :: z
complex             :: BesselY_RC

real    :: vm
complex :: CBJ, CDJ, CBY, CDY

call CJYVB(v, z, vm, CBJ, CDJ, CBY, CDY)

BesselY_RC = CBY

end function BesselY_RC

!-----------------------------------------------------------------------------!

elemental function BesselY_CR(zz,v)
real,    intent(in) :: v
complex, intent(in) :: zz
complex             :: BesselY_CR

complex :: BJ_1, BJ_2, z

z=CMPLX(v,0d0)

call cbsslj(zz,z, BJ_1)
call cbsslj(-zz,z, BJ_2)

BesselY_CR = (BJ_1*cos(zz*Pi)-BJ_2)/sin(zz*Pi)

end function BesselY_CR

!-----------------------------------------------------------------------------!

elemental function BesselY_CC(zz,z)
complex, intent(in) :: zz
complex, intent(in) :: z
complex             :: BesselY_CC

complex :: BJ_1, BJ_2

call cbsslj(zz,z, BJ_1)
call cbsslj(-zz,z, BJ_2)

BesselY_CC = (BJ_1*cos(zz*Pi)-BJ_2)/sin(zz*Pi)

end function BesselY_CC


!-----------------------------------------------------------------------------!
!                                   BesselI                                   !
!                                                                             !
! This set of functions returns the value of the modified Bessel function of  !
! the first kind I_n(x) for real as well as integer orders n and for real as  !
! well as complex arguments x.                                                !
!                                                                             !
! Inputs:                                                                     !
!    n: Order of the Bessel function.                                         !
!    x: Argument of the Bessel function.                                      !
!                                                                             !
! Outputs:                                                                    !
!    BesselI: Value of the Bessel function I_n(x).                            !
!                                                                             !
!-----------------------------------------------------------------------------!

elemental function BesselI_IR(n,x)
integer, intent(in) :: n
real,    intent(in) :: x
real                :: BesselI_IR

integer :: nm
real    :: BI, DI, BK, DK

call IKNB(n, x, nm, BI, DI, BK, DK)

BesselI_IR = BI

end function BesselI_IR

!-----------------------------------------------------------------------------!

elemental function BesselI_RR(v,x)
real,    intent(in) :: v
real,    intent(in) :: x
real                :: BesselI_RR

real    :: vm, BI, DI, BK, DK

call IKV(v, x, vm, BI, DI, BK, DK)

BesselI_RR = BI

end function BesselI_RR

!-----------------------------------------------------------------------------!

elemental function BesselI_IC(n,z)
integer, intent(in) :: n
complex, intent(in) :: z
complex             :: BesselI_IC

integer :: nm
complex :: CBI, CDI, CBK, CDK

call CIKNB(n, z, nm, CBI, CDI, CBK, CDK)

BesselI_IC = CBI

end function BesselI_IC

!-----------------------------------------------------------------------------!

elemental function BesselI_RC(v,z)
real,    intent(in) :: v
complex, intent(in) :: z
complex             :: BesselI_RC

real    :: vm
complex :: CBI, CDI, CBK, CDK

call CIKVB(v, z, vm, CBI, CDI, CBK, CDK)

BesselI_RC = CBI

end function BesselI_RC

!-----------------------------------------------------------------------------!

elemental function BesselI_CR(zz,v)
complex, intent(in) :: zz
real,    intent(in) :: v
complex             :: BesselI_CR

complex :: z, sol, i

i=CMPLX(0d0,1)

z=CMPLX(v,0d0)*i

call cbsslj(zz ,z, sol)

BesselI_CR = sol/(exp(zz*i*Pi/2))

end function BesselI_CR

!-----------------------------------------------------------------------------!

elemental function BesselI_CC(zz,z)
complex, intent(in) :: zz
complex, intent(in) :: z
complex             :: BesselI_CC

complex :: sol, i

i=CMPLX(0d0,1)

call cbsslj(zz ,z*i, sol)

BesselI_CC = sol/(exp(zz*i*Pi/2))

end function BesselI_CC


!-----------------------------------------------------------------------------!
!                                   BesselK                                   !
!                                                                             !
! This set of functions returns the value of the modified Bessel function of  !
! the second kind K_n(x) for real as well as integer orders n and for real as !
! well as complex arguments x.                                                !
!                                                                             !
! Inputs:                                                                     !
!    n: Order of the Bessel function.                                         !
!    x: Argument of the Bessel function.                                      !
!                                                                             !
! Outputs:                                                                    !
!    BesselK: Value of the Bessel function K_n(x).                            !
!                                                                             !
!-----------------------------------------------------------------------------!


elemental function BesselK_IR(n,x)
integer, intent(in) :: n
real,    intent(in) :: x
real                :: BesselK_IR

integer :: nm
real    :: BI, DI, BK, DK

call IKNB(n, x, nm, BI, DI, BK, DK)

BesselK_IR = BK

end function BesselK_IR

!-----------------------------------------------------------------------------!

elemental function BesselK_RR(v,x)
real,    intent(in) :: v
real,    intent(in) :: x
real                :: BesselK_RR

real    :: vm, BI, DI, BK, DK

call IKV(v, x, vm, BI, DI, BK, DK)

BesselK_RR = BK

end function BesselK_RR

!-----------------------------------------------------------------------------!

elemental function BesselK_IC(n,z)
integer, intent(in) :: n
complex, intent(in) :: z
complex             :: BesselK_IC

integer :: nm
complex :: CBI, CDI, CBK, CDK

call CIKNB(n, z, nm, CBI, CDI, CBK, CDK)

BesselK_IC = CBK

end function BesselK_IC

!-----------------------------------------------------------------------------!

elemental function BesselK_RC(v,z)
real,    intent(in) :: v
complex, intent(in) :: z
complex             :: BesselK_RC

real    :: vm
complex :: CBI, CDI, CBK, CDK

call CIKVB(v, z, vm, CBI, CDI, CBK, CDK)

BesselK_RC = CBK

end function BesselK_RC

!-----------------------------------------------------------------------------!

elemental function BesselK_CR(zz,v)
complex, intent(in) :: zz
real,    intent(in) :: v
complex             :: BesselK_CR

complex:: BI_1, BI_2

BI_1=BesselI(-zz, v)
BI_2=BesselI( zz, v)

BesselK_CR = (BI_1 -  BI_2)*(Pi/(2*sin(zz*Pi)))

end function BesselK_CR

!-----------------------------------------------------------------------------!

elemental function BesselK_CC(zz,z)
complex, intent(in) :: zz
complex, intent(in) :: z
complex             :: BesselK_CC

complex:: BI_1, BI_2

BI_1=BesselI(-zz, z)
BI_2=BesselI( zz, z)

BesselK_CC = (BI_1 -  BI_2)*(Pi/(2*sin(zz*Pi)))

end function BesselK_CC

!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!                                                                       !
!        Functions and subroutines that compute Bessel Functions        !
!                                                                       !
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------------!

elemental SUBROUTINE jynb(n, x, nm, bjj, djj, byy, dyy)

!    =====================================================
!    Purpose: Compute Bessel functions Jn(x), Yn(x) and
!             their derivatives
!    Input :  x --- Argument of Jn(x) and Yn(x) ( x Ú 0 )
!             n --- Order of Jn(x) and Yn(x)
!    Output:  BJ(n) --- Jn(x)
!             DJ(n) --- Jn'(x)
!             BY(n) --- Yn(x)
!             DY(n) --- Yn'(x)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 to calculate the starting
!             point for backward recurrence
!    =====================================================

INTEGER, INTENT(IN)     :: n
REAL,    INTENT(IN)     :: x
INTEGER, INTENT(OUT)    :: nm
REAL , INTENT(OUT)  :: bjj
REAL , INTENT(OUT)  :: djj
REAL , INTENT(OUT)  :: byy
REAL , INTENT(OUT)  :: dyy

REAL  :: bj(0:250)
REAL  :: dj(0:250)
REAL  :: by(0:250)
REAL  :: dy(0:250)

REAL , PARAMETER  ::  a(4) = (/ -.7031250000000000D-01,  &
     .1121520996093750D+00, -.5725014209747314D+00, .6074042001273483D+01 /)
REAL , PARAMETER  ::  b(4) = (/ .7324218750000000D-01,  &
    -.2271080017089844D+00, .1727727502584457D+01, -.2438052969955606D+02 /)
REAL , PARAMETER  ::  a1(4) = (/ .1171875000000000D+00,  &
    -.1441955566406250D+00, .6765925884246826D+00, -.6883914268109947D+01 /)
REAL , PARAMETER  ::  b1(4) = (/ -.1025390625000000D+00,  &
     .2775764465332031D+00, -.1993531733751297D+01, .2724882731126854D+02 /)
REAL , PARAMETER  :: pi = 3.141592653589793d0, r2p = .63661977236758d0
REAL   :: bj0, bj1, bjk, bs, by0, by1, byk, cu, ec, f, f1, f2,  &
              p0, p1, q0, q1, s0, su, sv, t1, t2
INTEGER    :: k, m

nm = n
IF (x < smallest_number) THEN  
  DO  k = 0, n
    bj(k) = 0.0D0
    dj(k) = 0.0D0
    by(k) = -HUGE(by(k))
    dy(k) = HUGE(dy(k))
  END DO
  bj(0) = 1.0D0
  dj(1) = 0.5D0
  RETURN
END IF
IF (x <= 300.0 .OR. n > INT(0.9*x)) THEN
  IF (n == 0) nm = 1
  m = msta1(x, 200)
  IF (m < nm) THEN
    nm = m
  ELSE
    m = msta2(x, nm, 15)
  END IF
  bs = 0.0D0
  su = 0.0D0
  sv = 0.0D0
  f2 = 0.0D0
  f1 = smallest_number
  DO  k = m, 0, -1
    f = 2.0D0 * (k+1) / x * f1 - f2
    IF (k <= nm) bj(k) = f
    IF (k == 2*INT(k/2) .AND. k /= 0) THEN
      bs = bs + 2.0D0 * f
      su = su + (-1) ** (k/2) * f / k
    ELSE IF (k > 1) THEN
      sv = sv + (-1) ** (k/2) * k / (k*k-1.0) * f
    END IF
    f2 = f1
    f1 = f
  END DO
  s0 = bs + f
  bj(0:nm) = bj(0:nm) / s0
  ec = LOG(x/2.0D0) + 0.5772156649015329d0
  by0 = r2p * (ec*bj(0)-4.0D0*su/s0)
  by(0) = by0
  by1 = r2p * ((ec-1.0D0)*bj(1) - bj(0)/x - 4.0D0*sv/s0)
  by(1) = by1
ELSE
  t1 = x - 0.25D0 * pi
  p0 = 1.0D0
  q0 = -0.125D0 / x
  DO  k = 1, 4
    p0 = p0 + a(k) * x ** (-2*k)
    q0 = q0 + b(k) * x ** (-2*k-1)
  END DO
  cu = SQRT(r2p/x)
  bj0 = cu * (p0*COS(t1) - q0*SIN(t1))
  by0 = cu * (p0*SIN(t1) + q0*COS(t1))
  bj(0) = bj0
  by(0) = by0
  t2 = x - 0.75D0 * pi
  p1 = 1.0D0
  q1 = 0.375D0 / x
  DO  k = 1, 4
    p1 = p1 + a1(k) * x ** (-2*k)
    q1 = q1 + b1(k) * x ** (-2*k-1)
  END DO
  bj1 = cu * (p1*COS(t2) - q1*SIN(t2))
  by1 = cu * (p1*SIN(t2) + q1*COS(t2))
  bj(1) = bj1
  by(1) = by1
  DO  k = 2, nm
    bjk = 2.0D0 * (k-1.0D0) / x * bj1 - bj0
    bj(k) = bjk
    bj0 = bj1
    bj1 = bjk
  END DO
END IF
dj(0) = -bj(1)
DO  k = 1, nm
  dj(k) = bj(k-1) - k / x * bj(k)
END DO
DO  k = 2, nm
  byk = 2.0D0 * (k-1.0D0) * by1 / x - by0
  by(k) = byk
  by0 = by1
  by1 = byk
END DO
dy(0) = -by(1)
DO  k = 1, nm
  dy(k) = by(k-1) - k * by(k) / x
END DO


bjj=bj(n)
djj=dj(n)
byy=by(n)
dyy=dy(n)

RETURN
END SUBROUTINE jynb

!-----------------------------------------------------------------------------!


elemental FUNCTION msta1(x, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that the magnitude of
!             Jn(x) at that point is about 10^(-MP)
!    Input :  x     --- Argument of Jn(x)
!             MP    --- Value of magnitude
!    Output:  MSTA1 --- Starting point
!    ===================================================

REAL , INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL   :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0d0 - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1

!-----------------------------------------------------------------------------!


elemental FUNCTION msta2(x, n, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that all Jn(x) has MP
!             significant digits
!    Input :  x  --- Argument of Jn(x)
!             n  --- Order of Jn(x)
!             MP --- Significant digit
!    Output:  MSTA2 --- Starting point
!    ===================================================

REAL , INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL   :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5d0 * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0d0 - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2

!-----------------------------------------------------------------------------!


elemental FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)    :: n
REAL , INTENT(IN)  :: x
REAL               :: fn_val

fn_val = 0.5d0 * LOG10(6.28d0*n) - n * LOG10(1.36d0*x/n)
RETURN
END FUNCTION envj

 
!-----------------------------------------------------------------------------!


elemental SUBROUTINE jyv(v, x, vm, bjj, djj, byy, dyy)

!    =======================================================
!    Purpose: Compute Bessel functions Jv(x) and Yv(x)
!             and their derivatives
!    Input :  x --- Argument of Jv(x) and Yv(x)
!             v --- Order of Jv(x) and Yv(x)
!                   ( v = n+v0, 0 Û v0 < 1, n = 0,1,2,... )
!    Output:  BJ(n) --- Jn+v0(x)
!             DJ(n) --- Jn+v0'(x)
!             BY(n) --- Yn+v0(x)
!             DY(n) --- Yn+v0'(x)
!             VM --- Highest order computed
!    Routines called:
!         (1) GAMMA for computing gamma function
!         (2) MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!    =======================================================

REAL, INTENT(IN)   :: v
REAL, INTENT(IN)   :: x
REAL, INTENT(OUT)  :: vm
REAL, INTENT(OUT)  :: bjj
REAL, INTENT(OUT)  :: djj
REAL, INTENT(OUT)  :: byy
REAL, INTENT(OUT)  :: dyy

REAL :: bj(0:250), dj(0:250), by(0:250), dy(0:250)

REAL, PARAMETER  :: el = .5772156649015329d0, pi = 3.141592653589793d0, &
                         rp2 = .63661977236758d0
REAL :: a, a0, b, bju0, bju1, bjv0, bjv1, bjvl, byv0, byv1, byvk,  &
              ck, cs, cs0, cs1, ec, f, f0, f1, f2, ga, gb, pv0, pv1, px,  &
              qx, r, r0, r1, rp, rq, sk, v0, vg, vl, vv, w0, w1, x2, xk
INTEGER    :: j, k, k0, l, m, n

x2 = x * x
n = v
v0 = v - n
IF (x < smallest_number) THEN  !FLAG
  DO  k = 0, n
    bj(k) = 0.0D0
    dj(k) = 0.0D0
    by(k) = -biggest_number
    dy(k) = biggest_number
  END DO
  IF (v0 == 0.0) THEN
    bj(0) = 1.0D0
    dj(1) = 0.5D0
  ELSE
    dj(0) = biggest_number
  END IF
  vm = v
  RETURN
END IF
IF (x <= 12.0) THEN
  DO  l = 0, 1
    vl = v0 + l
    bjvl = 1.0D0
    r = 1.0D0
    DO  k = 1, 40
      r = -0.25D0 * r * x2 / (k*(k+vl))
      bjvl = bjvl + r
      IF (ABS(r) < ABS(bjvl)*1.0D-15) EXIT
    END DO
    vg = 1.0D0 + vl
    CALL gamma(vg, ga)
    a = (0.5D0*x) ** vl / ga
    IF (l == 0) bjv0 = bjvl * a
    IF (l == 1) bjv1 = bjvl * a
  END DO
ELSE
  k0 = 11
  IF (x >= 35.0) k0 = 10
  IF (x >= 50.0) k0 = 8
  DO  j = 0, 1
    vv = 4.0D0 * (j+v0) * (j+v0)
    px = 1.0D0
    rp = 1.0D0
    DO  k = 1, k0
      rp = -0.78125D-2 * rp * (vv - (4*k-3)**2) * (vv - (4*k-1)**2) /  &
           (k*(2*k-1)*x2)
      px = px + rp
    END DO
    qx = 1.0D0
    rq = 1.0D0
    DO  k = 1, k0
      rq = -0.78125D-2 * rq * (vv - (4*k-1)**2) * (vv - (4*k+1)**2) /   &
           (k*(2*k+1)*x2)
      qx = qx + rq
    END DO
    qx = 0.125D0 * (vv-1.0) * qx / x
    xk = x - (0.5D0*(j+v0) + 0.25D0) * pi
    a0 = SQRT(rp2/x)
    ck = COS(xk)
    sk = SIN(xk)
    IF (j == 0) THEN
      bjv0 = a0 * (px*ck - qx*sk)
      byv0 = a0 * (px*sk + qx*ck)
    ELSE IF (j == 1) THEN
      bjv1 = a0 * (px*ck - qx*sk)
      byv1 = a0 * (px*sk + qx*ck)
    END IF
  END DO
END IF
bj(0) = bjv0
bj(1) = bjv1
dj(0) = v0 / x * bj(0) - bj(1)
dj(1) = -(1.0D0+v0) / x * bj(1) + bj(0)
IF (n >= 2 .AND. n <= INT(0.9*x)) THEN
  f0 = bjv0
  f1 = bjv1
  DO  k = 2, n
    f = 2.0D0 * (k+v0-1.0D0) / x * f1 - f0
    bj(k) = f
    f0 = f1
    f1 = f
  END DO
ELSE IF (n >= 2) THEN
  m = msta1(x, 200)
  IF (m < n) THEN
    n = m
  ELSE
    m = msta2(x, n, 15)
  END IF
  f2 = 0.0D0
  f1 = smallest_number
  DO  k = m, 0, -1
    f = 2.0D0 * (v0+k+1) / x * f1 - f2
    IF (k <= n) bj(k) = f
    f2 = f1
    f1 = f
  END DO
  IF (ABS(bjv0) > ABS(bjv1)) THEN
    cs = bjv0 / f
  ELSE
    cs = bjv1 / f2
  END IF
  bj(0:n) = cs * bj(0:n)
END IF
DO  k = 2, n
  dj(k) = -(k+v0) / x * bj(k) + bj(k-1)
END DO
IF (x <= 12.0D0) THEN
  IF (v0 /= 0.0) THEN
    DO  l = 0, 1
      vl = v0 + l
      bjvl = 1.0D0
      r = 1.0D0
      DO  k = 1, 40
        r = -0.25D0 * r * x2 / (k*(k-vl))
        bjvl = bjvl + r
        IF (ABS(r) < ABS(bjvl)*1.0D-15) EXIT
      END DO
      vg = 1.0D0 - vl
      CALL gamma(vg,gb)
      b = (2.0D0/x) ** vl / gb
      IF (l == 0) bju0 = bjvl * b
      IF (l == 1) bju1 = bjvl * b
    END DO
    pv0 = pi * v0
    pv1 = pi * (1.0D0+v0)
    byv0 = (bjv0*COS(pv0) - bju0) / SIN(pv0)
    byv1 = (bjv1*COS(pv1) - bju1) / SIN(pv1)
  ELSE
    ec = LOG(x/2.0D0) + el
    cs0 = 0.0D0
    w0 = 0.0D0
    r0 = 1.0D0
    DO  k = 1, 30
      w0 = w0 + 1.0D0 / k
      r0 = -0.25D0 * r0 / (k*k) * x2
      cs0 = cs0 + r0 * w0
    END DO
    byv0 = rp2 * (ec*bjv0 - cs0)
    cs1 = 1.0D0
    w1 = 0.0D0
    r1 = 1.0D0
    DO  k = 1, 30
      w1 = w1 + 1.0D0 / k
      r1 = -0.25D0 * r1 / (k*(k+1)) * x2
      cs1 = cs1 + r1 * (2.0D0*w1+1.0D0/(k+1.0D0))
    END DO
    byv1 = rp2 * (ec*bjv1 - 1.0D0/x - 0.25D0*x*cs1)
  END IF
END IF
by(0) = byv0
by(1) = byv1
DO  k = 2, n
  byvk = 2.0D0 * (v0+k-1) / x * byv1 - byv0
  by(k) = byvk
  byv0 = byv1
  byv1 = byvk
END DO
dy(0) = v0 / x * by(0) - by(1)
DO  k = 1, n
  dy(k) = -(k+v0) / x * by(k) + by(k-1)
END DO
vm = n + v0


bjj=bj(int(v))
djj=dj(int(v))
byy=by(int(v))
dyy=dy(int(v))


RETURN
END SUBROUTINE jyv

!-----------------------------------------------------------------------------!


elemental SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute gamma function ‚(x)
!    Input :  x  --- Argument of ‚(x)
!                    ( x is not equal to 0,-1,-2,˙˙˙)
!    Output:  GA --- ‚(x)
!    ==================================================

REAL, INTENT(IN)   :: x
REAL, INTENT(OUT)  :: ga

REAL, PARAMETER  :: g(26) = (/ 1.0d0, 0.5772156649015329d0,  &
     -0.6558780715202538d0, -0.420026350340952D-1, 0.1665386113822915d0,   &
     -0.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,  &
     -0.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3,  &
     -0.201348547807D-4, -.12504934821D-5, .11330272320D-5, -.2056338417D-6, &
      0.61160950D-8, .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /)
REAL, PARAMETER  :: pi = 3.141592653589793d0
REAL  :: gr, r, z
INTEGER    :: k, m, m1

IF (x == INT(x)) THEN
  IF (x > 0.0d0) THEN
    ga = 1.0d0
    m1 = x - 1
    DO  k = 2, m1
      ga = ga * k
    END DO
  ELSE
    ga = biggest_number
  END IF
ELSE
  IF (ABS(x) > 1.0d0) THEN
    z = ABS(x)
    m = z
    r = 1.0d0
    DO  k = 1, m
      r = r * (z-k)
    END DO
    z = z - m
  ELSE
    z = x
  END IF
  gr = g(26)
  DO  k = 25, 1, -1
    gr = gr * z + g(k)
  END DO
  ga = 1.0d0 / (gr*z)
  IF (ABS(x) > 1.0d0) THEN
    ga = ga * r
    IF (x < 0.0d0) ga = -pi / (x*ga*SIN(pi*x))
  END IF
END IF
END SUBROUTINE gamma
 

!-----------------------------------------------------------------------------!


elemental SUBROUTINE cjynb(n, z, nm, cbjj, cdjj, cbyy, cdyy)

!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =======================================================

INTEGER, INTENT(IN)   :: n
COMPLEX, INTENT(IN)   :: z
INTEGER, INTENT(OUT)  :: nm
COMPLEX, INTENT(OUT)  :: cbjj
COMPLEX, INTENT(OUT)  :: cdjj
COMPLEX, INTENT(OUT)  :: cbyy
COMPLEX, INTENT(OUT)  :: cdyy

COMPLEX :: cbj(0:250), cdj(0:250), cby(0:250), cdy(0:250)

REAL, PARAMETER  :: a(4) = (/ -.7031250000000000D-01, .1121520996093750d0,  &
      -0.5725014209747314d0, .6074042001273483D+01 /)
REAL, PARAMETER  :: b(4) = (/ .7324218750000000D-01, -.2271080017089844d0,  &
      .1727727502584457D+01, -.2438052969955606D+02 /)
REAL, PARAMETER  :: a1(4) = (/ .1171875000000000d0, -.1441955566406250d0,  &
      .6765925884246826d0, -.6883914268109947D+01 /)
REAL, PARAMETER  :: b1(4) = (/ -.1025390625000000d0, .2775764465332031d0,  &
      -0.1993531733751297D+01, .2724882731126854D+02 /)

REAL, PARAMETER  :: el = 0.5772156649015329D0, pi = 3.141592653589793D0, &
                         r2p = .63661977236758D0
COMPLEX  :: cbj0, cbj1, cbjk, cbs, cby0, cby1, ce, cf, cf1, cf2,  &
                 cp0, cp1, cq0, cq1, cs0, csu, csv, ct1, ct2, cu, cyy
REAL     :: a0, y0
INTEGER       :: k, m

y0 = ABS(AIMAG(z))
a0 = ABS(z)
nm = n
IF (a0 < smallest_number) THEN 
  DO  k = 0, n
    cbj(k) = (0.0D0,0.0D0)
    cdj(k) = (0.0D0,0.0D0)
    cby(k) = -(biggest_number,0.0D0) 
    cdy(k) = (biggest_number,0.0D0) 
  END DO
  cbj(0) = (1.0D0,0.0D0)
  cdj(1) = (0.5D0,0.0D0)
  RETURN
END IF
IF (a0 <= 300.d0 .OR. n > INT(0.25*a0)) THEN
  IF (n == 0) nm = 1
  m = msta1(a0, 200)
  IF (m < nm) THEN
    nm = m
  ELSE
    m = msta2(a0, nm, 15)
  END IF
  cbs = (0.0D0,0.0D0)
  csu = (0.0D0,0.0D0)
  csv = (0.0D0,0.0D0)
  cf2 = (0.0D0,0.0D0)
  cf1 = (smallest_number,0.0D0) !FLAG
  DO  k = m, 0, -1
    cf = 2.0D0 * (k+1.0D0) / z * cf1 - cf2
    IF (k <= nm) cbj(k) = cf
    IF (k == 2*INT(k/2) .AND. k /= 0) THEN
      IF (y0 <= 1.0D0) THEN
        cbs = cbs + 2.0D0 * cf
      ELSE
        cbs = cbs + (-1) ** (k/2) * 2.0D0 * cf
      END IF
      csu = csu + (-1) ** (k/2) * cf / k
    ELSE IF (k > 1) THEN
      csv = csv + (-1) ** (k/2) * k / (k*k-1.0D0) * cf
    END IF
    cf2 = cf1
    cf1 = cf
  END DO
  IF (y0 <= 1.0D0) THEN
    cs0 = cbs + cf
  ELSE
    cs0 = (cbs+cf) / COS(z)
  END IF
  DO  k = 0, nm
    cbj(k) = cbj(k) / cs0
  END DO
  ce = LOG(z/2.0D0) + el
  cby(0) = r2p * (ce*cbj(0) - 4.0D0*csu/cs0)
  cby(1) = r2p * (-cbj(0)/z + (ce-1.0D0)*cbj(1) - 4.0D0*csv/cs0)
ELSE
  ct1 = z - 0.25D0 * pi
  cp0 = (1.0D0,0.0D0)
  DO  k = 1, 4
    cp0 = cp0 + a(k) * z ** (-2*k)
  END DO
  cq0 = -0.125D0 / z
  DO  k = 1, 4
    cq0 = cq0 + b(k) * z ** (-2*k-1)
  END DO
  cu = SQRT(r2p/z)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  cbj(0) = cbj0
  cby(0) = cby0
  ct2 = z - 0.75D0 * pi
  cp1 = (1.0D0,0.0D0)
  DO  k = 1, 4
    cp1 = cp1 + a1(k) * z ** (-2*k)
  END DO
  cq1 = 0.375D0 / z
  DO  k = 1, 4
    cq1 = cq1 + b1(k) * z ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
  cbj(1) = cbj1
  cby(1) = cby1
  DO  k = 2, nm
    cbjk = 2.0D0 * (k-1.0D0) / z * cbj1 - cbj0
    cbj(k) = cbjk
    cbj0 = cbj1
    cbj1 = cbjk
  END DO
END IF
cdj(0) = -cbj(1)
DO  k = 1, nm
  cdj(k) = cbj(k-1) - k / z * cbj(k)
END DO
IF (ABS(cbj(0)) > 1.0D0) THEN
  cby(1) = (cbj(1)*cby(0)-2.0D0/(pi*z)) / cbj(0)
END IF
DO  k = 2, nm
  IF (ABS(cbj(k-1)) >= ABS(cbj(k-2))) THEN
    cyy = (cbj(k)*cby(k-1)-2.0D0/(pi*z)) / cbj(k-1)
  ELSE
    cyy = (cbj(k)*cby(k-2)-4.0D0*(k-1.0D0)/(pi*z*z)) / cbj(k-2)
  END IF
  cby(k) = cyy
END DO
cdy(0) = -cby(1)
DO  k = 1, nm
  cdy(k) = cby(k-1) - k / z * cby(k)
END DO

cbjj=cbj(n)
cdjj=cdj(n)
cbyy=cby(n)
cdyy=cdy(n)


RETURN
END SUBROUTINE cjynb

!-----------------------------------------------------------------------------!


elemental SUBROUTINE cjyvb(v, z, vm, cbjj, cdjj, cbyy, cdyy)

!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 Û v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================

REAL, INTENT(IN)      :: v
COMPLEX, INTENT(IN)   :: z
REAL, INTENT(OUT)     :: vm
COMPLEX, INTENT(OUT)  :: cbjj
COMPLEX, INTENT(OUT)  :: cdjj
COMPLEX, INTENT(OUT)  :: cbyy
COMPLEX, INTENT(OUT)  :: cdyy

COMPLEX :: cbj(0:250), cdj(0:250), cby(0:250), cdy(0:250)

REAL, PARAMETER  :: pi = 3.141592653589793D0, rp2 = .63661977236758D0
COMPLEX  :: ca, ca0, cb, cck, cec, cf, cf1, cf2, cfac0, ci,  &
                 cju0, cjv0, cjvn, cpz, cqz, cr, cr0, crp, crq,  &
                 cs, cs0, csk, cyv0, cyy, z1, z2, zk
REAL     :: a0, ga, gb, pv0, v0, vg, vv, w0
INTEGER       :: k, k0, m, n

ci = (0.0D0,1.0D0)
a0 = ABS(z)
z1 = z
z2 = z * z
n = INT(v)
v0 = v - n
pv0 = pi * v0
IF (a0 < smallest_number) THEN 
  DO  k = 0, n
    cbj(k) = (0.0D0,0.0D0)
    cdj(k) = (0.0D0,0.0D0)
    cby(k) = -(biggest_number,0.0D0) 
    cdy(k) = (biggest_number,0.0D0) 
  END DO
  IF (v0 == 0.0) THEN
    cbj(0) = (1.0D0,0.0D0)
    cdj(1) = (0.5D0,0.0D0)
  ELSE
    cdj(0) = (biggest_number,0.0D0) 
  END IF
  vm = v
  RETURN
END IF
IF (REAL(z) < 0.0D0) z1 = -z
IF (a0 <= 12.0) THEN
  cjv0 = (1.0D0,0.0D0)
  cr = (1.0D0,0.0D0)
  DO  k = 1, 40
    cr = -0.25D0 * cr * z2 / (k*(k+v0))
    cjv0 = cjv0 + cr
    IF (ABS(cr) < ABS(cjv0)*1.0D-15) EXIT
  END DO

  vg = 1.0D0 + v0
  CALL gamma(vg, ga)
  ca = (0.5D0*z1) ** v0 / ga
  cjv0 = cjv0 * ca
ELSE
  k0 = 11
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  vv = 4.0D0 * v0 * v0
  cpz = (1.0D0,0.0D0)
  crp = (1.0D0,0.0D0)
  DO  k = 1, k0
    crp = -0.78125D-2 * crp * (vv - (4*k-3)**2) * (vv - (4*k-1)**2) / (k*(2*k-1)*z2)
    cpz = cpz + crp
  END DO
  cqz = (1.0D0,0.0D0)
  crq = (1.0D0,0.0D0)
  DO  k = 1, k0
    crq = -0.78125D-2 * crq * (vv - (4*k-1)**2) * (vv - (4*k+1)**2) / (k*(2*k+1)*z2)
    cqz = cqz + crq
  END DO
  cqz = 0.125D0 * (vv-1.0) * cqz / z1
  zk = z1 - (0.5D0*v0 + 0.25D0) * pi
  ca0 = SQRT(rp2/z1)
  cck = COS(zk)
  csk = SIN(zk)
  cjv0 = ca0 * (cpz*cck - cqz*csk)
  cyv0 = ca0 * (cpz*csk + cqz*cck)
END IF
IF (a0 <= 12.0) THEN
  IF (v0 /= 0.0) THEN
    cjvn = (1.0D0,0.0D0)
    cr = (1.0D0,0.0D0)
    DO  k = 1, 40
      cr = -0.25D0 * cr * z2 / (k*(k-v0))
      cjvn = cjvn + cr
      IF (ABS(cr) < ABS(cjvn)*1.0D-15) EXIT
    END DO

    vg = 1.0D0 - v0
    CALL gamma(vg,gb)
    cb = (2.0D0/z1) ** v0 / gb
    cju0 = cjvn * cb
    cyv0 = (cjv0*COS(pv0) - cju0) / SIN(pv0)
  ELSE
    cec = LOG(z1/2.0D0) + .5772156649015329D0
    cs0 = (0.0D0,0.0D0)
    w0 = 0.0D0
    cr0 = (1.0D0,0.0D0)
    DO  k = 1, 30
      w0 = w0 + 1.0D0 / k
      cr0 = -0.25D0 * cr0 / (k*k) * z2
      cs0 = cs0 + cr0 * w0
    END DO
    cyv0 = rp2 * (cec*cjv0-cs0)
  END IF
END IF
IF (n == 0) n = 1
m = msta1(a0, 200)
IF (m < n) THEN
  n = m
ELSE
  m = msta2(a0, n, 15)
END IF
cf2 = (0.0D0,0.0D0)
cf1 = (smallest_number,0.0D0)  
DO  k = m, 0, -1
  cf = 2.0D0 * (v0+k+1) / z1 * cf1 - cf2
  IF (k <= n) cbj(k) = cf
  cf2 = cf1
  cf1 = cf
END DO
cs = cjv0 / cf
DO  k = 0, n
  cbj(k) = cs * cbj(k)
END DO
IF (REAL(z) < 0.0D0) THEN
  cfac0 = EXP(pv0*ci)
  IF (AIMAG(z) < 0.0D0) THEN
    cyv0 = cfac0 * cyv0 - 2.0D0 * ci * COS(pv0) * cjv0
  ELSE IF (AIMAG(z) > 0.0D0) THEN
    cyv0 = cyv0 / cfac0 + 2.0D0 * ci * COS(pv0) * cjv0
  END IF
  DO  k = 0, n
    IF (AIMAG(z) < 0.0D0) THEN
      cbj(k) = EXP(-pi*(k+v0)*ci) * cbj(k)
    ELSE IF (AIMAG(z) > 0.0D0) THEN
      cbj(k) = EXP(pi*(k+v0)*ci) * cbj(k)
    END IF
  END DO
END IF
cby(0) = cyv0
DO  k = 1, n
  cyy = (cbj(k)*cby(k-1) - 2.0D0/(pi*z)) / cbj(k-1)
  cby(k) = cyy
END DO
cdj(0) = v0 / z * cbj(0) - cbj(1)
DO  k = 1, n
  cdj(k) = -(k+v0) / z * cbj(k) + cbj(k-1)
END DO
cdy(0) = v0 / z * cby(0) - cby(1)
DO  k = 1, n
  cdy(k) = cby(k-1) - (k+v0) / z * cby(k)
END DO
vm = n + v0


cbjj=cbj(int(v))
cdjj=cdj(int(v))
cbyy=cby(int(v))
cdyy=cdy(int(v))


RETURN
END SUBROUTINE cjyvb

!-----------------------------------------------------------------------------!

elemental SUBROUTINE iknb(n,x,nm,bii,dii,bkk,dkk)

!    ============================================================
!    Purpose: Compute modified Bessel functions In(x) and Kn(x),
!             and their derivatives
!    Input:   x --- Argument of In(x) and Kn(x) ( 0 Û x Û 700 )
!             n --- Order of In(x) and Kn(x)
!    Output:  BI(n) --- In(x)
!             DI(n) --- In'(x)
!             BK(n) --- Kn(x)
!             DK(n) --- Kn'(x)
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 for computing the starting point
!             for backward recurrence
!    ===========================================================

INTEGER, INTENT(IN)     :: n
REAL, INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL, INTENT(OUT)  :: bii
REAL, INTENT(OUT)  :: dii
REAL, INTENT(OUT)  :: bkk
REAL, INTENT(OUT)  :: dkk

REAL :: bi(0:250), di(0:250), bk(0:250), dk(0:250)

REAL, PARAMETER  :: pi = 3.141592653589793d0, el = 0.5772156649015329d0
REAL  :: a0, bkl, bs, f, f0, f1, g, g0, g1, r, s0, sk0, vt
INTEGER    :: k, k0, l, m

nm = n
IF (x <= smallest_number) THEN 
  DO  k = 0, n
    bi(k) = 0.0D0
    di(k) = 0.0D0
    bk(k) = biggest_number
    dk(k) = -biggest_number
  END DO
  bi(0) = 1.0D0
  di(1) = 0.5D0
  RETURN
END IF
IF (n == 0) nm = 1
m = msta1(x, 200)
IF (m < nm) THEN
  nm = m
ELSE
  m = msta2(x, nm, 15)
END IF
bs = 0.0D0
sk0 = 0.0D0
f0 = 0.0D0
f1 = smallest_number
DO  k = m, 0, -1
  f = 2*(k+1)/x*f1 + f0
  IF (k <= nm) bi(k) = f
  IF (k /= 0 .AND. k == 2*INT(k/2)) sk0 = sk0 + 4.0D0 * f / k
  bs = bs + 2.0D0 * f
  f0 = f1
  f1 = f
END DO
s0 = EXP(x) / (bs-f)
bi(0:nm) = s0 * bi(0:nm)
IF (x <= 8.0D0) THEN
  bk(0) = -(LOG(0.5D0*x)+el) * bi(0) + s0 * sk0
  bk(1) = (1.0D0/x-bi(1)*bk(0)) / bi(0)
ELSE
  a0 = SQRT(pi/(2.0D0*x)) * EXP(-x)
  k0 = 16
  IF (x >= 25.0) k0 = 10
  IF (x >= 80.0) k0 = 8
  IF (x >= 200.0) k0 = 6
  DO  l = 0, 1
    bkl = 1.0D0
    vt = 4 * l
    r = 1.0D0
    DO  k = 1, k0
      r = 0.125D0 * r * (vt - (2*k-1)**2) / (k*x)
      bkl = bkl + r
    END DO
    bk(l) = a0 * bkl
  END DO
END IF
g0 = bk(0)
g1 = bk(1)
DO  k = 2, nm
  g = 2*(k-1)/x*g1 + g0
  bk(k) = g
  g0 = g1
  g1 = g
END DO
di(0) = bi(1)
dk(0) = -bk(1)
DO  k = 1, nm
  di(k) = bi(k-1) - k / x * bi(k)
  dk(k) = -bk(k-1) - k / x * bk(k)
END DO

bii=bi(n)
dii=di(n)
bkk=bk(n)
dkk=dk(n)

RETURN

END SUBROUTINE iknb

!-----------------------------------------------------------------------------!


elemental SUBROUTINE ikv(v, x, vm, bii, dii, bkk, dkk)

!    =======================================================
!    Purpose: Compute modified Bessel functions Iv(x) and
!             Kv(x), and their derivatives
!    Input :  x --- Argument ( x Ú 0 )
!             v --- Order of Iv(x) and Kv(x)
!                   ( v = n+v0, n = 0,1,2,..., 0 Û v0 < 1 )
!    Output:  BI(n) --- In+v0(x)
!             DI(n) --- In+v0'(x)
!             BK(n) --- Kn+v0(x)
!             DK(n) --- Kn+v0'(x)
!             VM --- Highest order computed
!    Routines called:
!         (1) GAMMA for computing the gamma function
!         (2) MSTA1 and MSTA2 to compute the starting
!             point for backward recurrence
!    =======================================================

REAL, INTENT(IN)   :: v
REAL, INTENT(IN)   :: x
REAL, INTENT(OUT)  :: vm
REAL, INTENT(OUT)  :: bii
REAL, INTENT(OUT)  :: dii
REAL, INTENT(OUT)  :: bkk
REAL, INTENT(OUT)  :: dkk

real    :: BI(0:250), DI(0:250), BK(0:250), DK(0:250)

REAL, PARAMETER  :: pi = 3.141592653589793d0
REAL  :: a1, a2, bi0, bk0, bk1, bk2, ca, cb, cs, ct, f, f1, f2, gan, gap,  &
              piv, r, r1, r2, sum, v0, v0n, v0p, vt, w0, wa, ww, x2
INTEGER    :: k, k0, m, n

x2 = x * x
n = INT(v)
v0 = v - n
IF (n == 0) n = 1
IF (x < smallest_number) THEN 
  DO  k = 0, n
    bi(k) = 0.0D0
    di(k) = 0.0D0
    bk(k) = -biggest_number
    dk(k) = biggest_number
  END DO
  IF (v == 0.0) THEN
    bi(0) = 1.0D0
    di(1) = 0.5D0
  END IF
  vm = v
  RETURN
END IF
piv = pi * v0
vt = 4.0D0 * v0 * v0
IF (v0 == 0.0D0) THEN
  a1 = 1.0D0
ELSE
  v0p = 1.0D0 + v0
  CALL gamma(v0p, gap)
  a1 = (0.5D0*x) ** v0 / gap
END IF
k0 = 14
IF (x >= 35.0) k0 = 10
IF (x >= 50.0) k0 = 8
IF (x <= 18.0) THEN
  bi0 = 1.0D0
  r = 1.0D0
  DO  k = 1, 30
    r = 0.25D0 * r * x2 / (k*(k+v0))
    bi0 = bi0 + r
    IF (ABS(r/bi0) < 1.0D-15) EXIT
  END DO
  bi0 = bi0 * a1
ELSE
  ca = EXP(x) / SQRT(2.0D0*pi*x)
  sum = 1.0D0
  r = 1.0D0
  DO  k = 1, k0
    r = -0.125D0 * r * (vt - (2*k-1)**2) / (k*x)
    sum = sum + r
  END DO
  bi0 = ca * sum
END IF
m = msta1(x,200)
IF (m < n) THEN
  n = m
ELSE
  m = msta2(x, n, 15)
END IF
f2 = 0.0D0
f1 = smallest_number
DO  k = m, 0, -1
  f = 2.0D0 * (v0+k+1) / x * f1 + f2
  IF (k <= n) bi(k) = f
  f2 = f1
  f1 = f
END DO
cs = bi0 / f
DO  k = 0, n
  bi(k) = cs * bi(k)
END DO
di(0) = v0 / x * bi(0) + bi(1)
DO  k = 1, n
  di(k) = -(k+v0) / x * bi(k) + bi(k-1)
END DO
IF (x <= 9.0D0) THEN
  IF (v0 == 0.0D0) THEN
    ct = -LOG(0.5D0*x) - 0.5772156649015329d0
    cs = 0.0D0
    w0 = 0.0D0
    r = 1.0D0
    DO  k = 1, 50
      w0 = w0 + 1.0D0 / k
      r = 0.25D0 * r / (k*k) * x2
      cs = cs + r * (w0+ct)
      wa = ABS(cs)
      IF (ABS((wa-ww)/wa) < 1.0D-15) EXIT
      ww = wa
    END DO
    bk0 = ct + cs
  ELSE
    v0n = 1.0D0 - v0
    CALL gamma(v0n,gan)
    a2 = 1.0D0 / (gan*(0.5D0*x)**v0)
    a1 = (0.5D0*x) ** v0 / gap
    sum = a2 - a1
    r1 = 1.0D0
    r2 = 1.0D0
    DO  k = 1, 120
      r1 = 0.25D0 * r1 * x2 / (k*(k-v0))
      r2 = 0.25D0 * r2 * x2 / (k*(k+v0))
      sum = sum + a2 * r1 - a1 * r2
      wa = ABS(sum)
      IF (ABS((wa-ww)/wa) < 1.0D-15) EXIT
      ww = wa
    END DO
    bk0 = 0.5D0 * pi * sum / SIN(piv)
  END IF
ELSE
  cb = EXP(-x) * SQRT(0.5D0*pi/x)
  sum = 1.0D0
  r = 1.0D0
  DO  k = 1, k0
    r = 0.125D0 * r * (vt-(2.0*k-1.0)**2.0) / (k*x)
    sum = sum + r
  END DO
  bk0 = cb * sum
END IF
bk1 = (1.0D0/x-bi(1)*bk0) / bi(0)
bk(0) = bk0
bk(1) = bk1
DO  k = 2, n
  bk2 = 2.0D0 * (v0+k-1.0D0) / x * bk1 + bk0
  bk(k) = bk2
  bk0 = bk1
  bk1 = bk2
END DO
dk(0) = v0 / x * bk(0) - bk(1)
DO  k = 1, n
  dk(k) = -(k+v0) / x * bk(k) - bk(k-1)
END DO
vm = n + v0

bii=bi(int(v))
dii=di(int(v))
bkk=bk(int(v))
dkk=dk(int(v))

RETURN

END SUBROUTINE ikv

!-----------------------------------------------------------------------------!


elemental SUBROUTINE ciknb(n, z, nm, cbii, cdii, cbkk, cdkk)

!       ============================================================
!       Purpose: Compute modified Bessel functions In(z) and Kn(z),
!                and their derivatives for a complex argument
!       Input:   z --- Complex argument
!                n --- Order of In(z) and Kn(z)
!       Output:  CBI(n) --- In(z)
!                CDI(n) --- In'(z)
!                CBK(n) --- Kn(z)
!                CDK(n) --- Kn'(z)
!                NM --- Highest order computed
!       Routones called:
!                MSTA1 and MSTA2 to compute the starting point for
!                backward recurrence
!       ===========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX, INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX, INTENT(OUT)  :: cbii
COMPLEX, INTENT(OUT)  :: cdii
COMPLEX, INTENT(OUT)  :: cbkk
COMPLEX, INTENT(OUT)  :: cdkk

COMPLEX:: cbi(0:250), cdi(0:250), cbk(0:250), cdk(0:250)

REAL, PARAMETER  :: pi = 3.141592653589793d0, el = 0.57721566490153d0
REAL     :: a0, fac, vt
INTEGER       :: k, k0, l, m
COMPLEX  :: ca0, cbkl, cbs, cf, cf0, cf1, cg, cg0, cg1, ci, cr, cs0, csk0, z1

a0 = ABS(z)
nm = n
IF (a0 < smallest_number) THEN  
  DO  k = 0, n
    cbi(k) = (0.0d0,0.0d0)
    cbk(k) = (biggest_number,0.0d0) 
    cdi(k) = (biggest_number,0.0d0)
    cdk(k) = -(biggest_number,0.0d0) 
  END DO
  cbi(0) = (1.0d0,0.0d0)
  cdi(1) = (0.5d0,0.0d0)
  RETURN
END IF
z1 = z
ci = (0.0d0,1.0d0)
IF (REAL(z) < 0.0) z1 = -z
IF (n == 0) nm = 1
m = msta1(a0, 200)
IF (m < nm) THEN
  nm = m
ELSE
  m = msta2(a0, nm, 15)
END IF
cbs = 0.0d0
csk0 = 0.0d0
cf0 = 0.0d0
cf1 = smallest_number
DO  k = m, 0, -1
  cf = 2.0d0 * (k+1.0d0) * cf1 / z1 + cf0
  IF (k <= nm) cbi(k) = cf
  IF (k /= 0 .AND. k == 2*INT(k/2)) csk0 = csk0 + 4.0d0 * cf / k
  cbs = cbs + 2.0d0 * cf
  cf0 = cf1
  cf1 = cf
END DO
cs0 = EXP(z1) / (cbs-cf)
DO  k = 0, nm
  cbi(k) = cs0 * cbi(k)
END DO
IF (a0 <= 9.0) THEN
  cbk(0) = -(LOG(0.5d0*z1)+el) * cbi(0) + cs0 * csk0
  cbk(1) = (1.0d0/z1-cbi(1)*cbk(0)) / cbi(0)
ELSE
  ca0 = SQRT(pi/(2.0d0*z1)) * EXP(-z1)
  k0 = 16
  IF (a0 >= 25.0) k0 = 10
  IF (a0 >= 80.0) k0 = 8
  IF (a0 >= 200.0) k0 = 6
  DO  l = 0, 1
    cbkl = 1.0d0
    vt = 4.0d0 * l
    cr = (1.0d0,0.0d0)
    DO  k = 1, k0
      cr = 0.125d0 * cr * (vt-(2.0*k-1.0)**2) / (k*z1)
      cbkl = cbkl + cr
    END DO
    cbk(l) = ca0 * cbkl
  END DO
END IF
cg0 = cbk(0)
cg1 = cbk(1)
DO  k = 2, nm
  cg = 2.0d0 * (k-1.0d0) / z1 * cg1 + cg0
  cbk(k) = cg
  cg0 = cg1
  cg1 = cg
END DO
IF (REAL(z) < 0.0) THEN
  fac = 1.0d0
  DO  k = 0, nm
    IF (AIMAG(z) < 0.0) THEN
      cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
    ELSE
      cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
    END IF
    cbi(k) = fac * cbi(k)
    fac = -fac
  END DO
END IF
cdi(0) = cbi(1)
cdk(0) = -cbk(1)
DO  k = 1, nm
  cdi(k) = cbi(k-1) - k / z * cbi(k)
  cdk(k) = -cbk(k-1) - k / z * cbk(k)
END DO

cbii=cbi(n)
cdii=cdi(n)
cbkk=cbk(n)
cdkk=cdk(n)

RETURN
END SUBROUTINE ciknb

!-----------------------------------------------------------------------------!


elemental SUBROUTINE cikvb(v, z, vm, cbii, cdii, cbkk, cdkk)

!       ===========================================================
!       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
!                and their derivatives for an arbitrary order and
!                complex argument
!       Input :  z --- Complex argument z
!                v --- Real order of Iv(z) and Kv(z)
!                      ( v =n+v0, n = 0,1,2,..., 0 Û v0 < 1 )
!       Output:  CBI(n) --- In+v0(z)
!                CDI(n) --- In+v0'(z)
!                CBK(n) --- Kn+v0(z)
!                CDK(n) --- Kn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================

REAL, INTENT(IN)      :: v
COMPLEX, INTENT(IN)   :: z
REAL, INTENT(OUT)     :: vm
COMPLEX, INTENT(OUT)  :: cbii
COMPLEX, INTENT(OUT)  :: cdii
COMPLEX, INTENT(OUT)  :: cbkk
COMPLEX, INTENT(OUT)  :: cdkk

COMPLEX:: cbi(0:250), cdi(0:250), cbk(0:250), cdk(0:250)

COMPLEX  :: ca, ca1, ca2, cb, cbi0, cbk0, cf, cf1, cf2, ci, ci0, ckk,  &
                 cp, cr, cr1, cr2, cs, csu, ct, cvk, z1, z2
REAL     :: a0, gan, gap, piv, v0, v0n, v0p, vt, w0
INTEGER       :: k, k0, m, n
REAL, PARAMETER  :: pi = 3.141592653589793d0

z1 = z
z2 = z * z
a0 = ABS(z)
ci = (0.0d0,1.0d0)
n = INT(v)
v0 = v - n
piv = pi * v0
vt = 4.0d0 * v0 * v0
IF (n == 0) n = 1
IF (a0 < 1.0) THEN 
  DO  k = 0, n
    cbi(k) = 0.0d0
    cdi(k) = 0.0d0
    cbk(k) = -biggest_number
    cdk(k) = biggest_number
  END DO
  IF (v0 == 0.0) THEN
    cbi(0) = (1.0d0,0.0d0)
    cdi(1) = (0.5d0,0.0d0)
  END IF
  vm = v
  RETURN
END IF
k0 = 14
IF (a0 >= 35.0) k0 = 10
IF (a0 >= 50.0) k0 = 8
IF (REAL(z) < 0.0) z1 = -z
IF (a0 < 18.0) THEN
  IF (v0 == 0.0) THEN
    ca1 = (1.0d0,0.0d0)
  ELSE
    v0p = 1.0d0 + v0
    CALL gamma(v0p, gap)
    ca1 = (0.5d0*z1) ** v0 / gap
  END IF
  ci0 = (1.0d0,0.0d0)
  cr = (1.0d0,0.0d0)
  DO  k = 1, 50
    cr = 0.25d0 * cr * z2 / (k*(k+v0))
    ci0 = ci0 + cr
    IF (ABS(cr/ci0) < 1.0D-15) EXIT
  END DO
  cbi0 = ci0 * ca1
ELSE
  ca = EXP(z1) / SQRT(2.0d0*pi*z1)
  cs = (1.0d0,0.0d0)
  cr = (1.0d0,0.0d0)
  DO  k = 1, k0
    cr = -0.125d0 * cr * (vt-(2.0d0*k-1.0d0)**2.0) / (k*z1)
    cs = cs + cr
  END DO
  cbi0 = ca * cs
END IF
m = msta1(a0, 200)
IF (m < n) THEN
  n = m
ELSE
  m = msta2(a0, n, 15)
END IF
cf2 = (0.0d0,0.0d0)
cf1 = (smallest_number,0.0d0)
DO  k = m, 0, -1
  cf = 2.0d0 * (v0+k+1) / z1 * cf1 + cf2
  IF (k <= n) cbi(k) = cf
  cf2 = cf1
  cf1 = cf
END DO
cs = cbi0 / cf
cbi(0:n) = cs * cbi(0:n)
IF (a0 <= 9.0) THEN
  IF (v0 == 0.0) THEN
    ct = -LOG(0.5d0*z1) - 0.5772156649015329d0
    cs = (0.0d0,0.0d0)
    w0 = 0.0d0
    cr = (1.0d0,0.0d0)
    DO  k = 1, 50
      w0 = w0 + 1.0d0 / k
      cr = 0.25d0 * cr / (k*k) * z2
      cp = cr * (w0+ct)
      cs = cs + cp
      IF (k >= 10 .AND. ABS(cp/cs) < 1.0D-15) EXIT
    END DO
    cbk0 = ct + cs
  ELSE
    v0n = 1.0d0 - v0
    CALL gamma(v0n,gan)
    ca2 = 1.0d0 / (gan*(0.5d0*z1)**v0)
    ca1 = (0.5d0*z1) ** v0 / gap
    csu = ca2 - ca1
    cr1 = (1.0d0,0.0d0)
    cr2 = (1.0d0,0.0d0)
    DO  k = 1, 50
      cr1 = 0.25d0 * cr1 * z2 / (k*(k-v0))
      cr2 = 0.25d0 * cr2 * z2 / (k*(k+v0))
      cp = ca2 * cr1 - ca1 * cr2
      csu = csu + cp
      IF (k >= 10 .AND. ABS(cp/csu) < 1.0D-15) EXIT
    END DO
    cbk0 = 0.5d0 * pi * csu / SIN(piv)
  END IF
ELSE
  cb = EXP(-z1) * SQRT(0.5d0*pi/z1)
  cs = (1.0d0,0.0d0)
  cr = (1.0d0,0.0d0)
  DO  k = 1, k0
    cr = 0.125d0 * cr * (vt-(2*k-1)**2.0) / (k*z1)
    cs = cs + cr
  END DO
  cbk0 = cb * cs
END IF
cbk(0) = cbk0
IF (REAL(z) < 0.0) THEN
  DO  k = 0, n
    cvk = EXP((k+v0)*pi*ci)
    IF (AIMAG(z) < 0.0d0) THEN
      cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
      cbi(k) = cbi(k) / cvk
    ELSE IF (AIMAG(z) > 0.0) THEN
      cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
      cbi(k) = cvk * cbi(k)
    END IF
  END DO
END IF
DO  k = 1, n
  ckk = (1.0d0/z - cbi(k)*cbk(k-1)) / cbi(k-1)
  cbk(k) = ckk
END DO
cdi(0) = v0 / z * cbi(0) + cbi(1)
cdk(0) = v0 / z * cbk(0) - cbk(1)
DO  k = 1, n
  cdi(k) = -(k+v0) / z * cbi(k) + cbi(k-1)
  cdk(k) = -(k+v0) / z * cbk(k) - cbk(k-1)
END DO
vm = n + v0

cbii=cbi(int(v))
cdii=cdi(int(v))
cbkk=cbk(int(v))
cdkk=cdk(int(v))


RETURN
END SUBROUTINE cikvb

!-----------------------------------------------------------------------------!



!-----------------------------------------------------------------------------!
!
!       COMPLEX FUNCTIONS FOR BESSELJ(C,C) FLAG
!
!-----------------------------------------------------------------------------!



elemental SUBROUTINE cbsslj(cnu,z,w)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEX BESSEL FUNCTION J   (Z)
!                                                    CNU
!-----------------------------------------------------------------------

!     WRITTEN BY
!         ANDREW H. VAN TUYL AND ALFRED H. MORRIS, JR.
!         NAVAL SURFACE WARFARE CENTER
!         OCTOBER, 1991

!     A MODIFICATION OF THE PROCEDURE DEVELOPED BY ALLEN V. HERSHEY
!     (NAVAL SURFACE WARFARE CENTER) IN 1978 FOR HANDLING THE DEBYE
!     APPROXIMATION IS EMPLOYED.

!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)   :: z
COMPLEX , INTENT(IN)   :: cnu
COMPLEX , INTENT(OUT)  :: w

COMPLEX   :: c, nu, s, sm1, sm2, t, tsc, w0, w1, zn, zz
!-----------------------
REAL  :: a, cn1, cn2, e, fn
REAL  :: pn, qm, qn, qnp1
REAL  :: r, rn2, r2, sn, t1, t2
REAL  :: u, v, x, y
INTEGER   :: i, k, m, n
REAL , PARAMETER  :: pi = 3.141592653589793238462643383279502884197d0
!-----------------------
x = REAL(z)
y = AIMAG(z)
r = cpabs(x,y)
cn1 = REAL(cnu)
cn2 = AIMAG(cnu)
rn2 = cn1 * cn1 + cn2 * cn2
pn = INT(cn1)
fn = cn1 - pn
sn = 1.0d0

!          CALCULATION WHEN ORDER IS AN INTEGER

IF (fn == 0.0d0 .AND. cn2 == 0.0d0) THEN
  n = pn
  pn = ABS(pn)
  cn1 = pn
  IF (n < 0 .AND. n /= (n/2)*2) sn = -1.0d0
END IF

!          SELECTION OF METHOD

IF (r > 17.5D0) THEN
  IF (r > 17.5D0 + 0.5D0*rn2) GO TO 10
  GO TO 20
END IF

!          USE MACLAURIN EXPANSION AND RECURSION

IF (cn1 < 0.0D0) THEN
  qn = -1.25D0 * (r + 0.5D0*ABS(cn2) - ABS(y-0.5D0*cn2))
  IF (cn1 < qn) THEN
    qn = 1.25D0 * (r - MAX(1.2D0*r,ABS(y-cn2)))
    IF (cn1 < qn) THEN
      qn = MIN(pn, REAL(-INT(1.25D0*(r-ABS(cn2)))))
      GO TO 60
    END IF
  END IF
END IF

r2 = r * r
qm = 0.0625D0 * r2 * r2 - cn2 * cn2
qn = MAX(pn, REAL(INT(SQRT(MAX(0.0D0,qm)))))
GO TO 60

!          USE ASYMPTOTIC EXPANSION

10 CALL cbja(z,cnu,w)
RETURN

!          CALCULATION FOR 17.5 < ABS(Z) <= 17.5 + 0.5*ABS(CNU)**2

20 n = 0
IF (ABS(cn2) < 0.8D0*ABS(y)) THEN
  qm = -1.25D0 * (r + 0.5D0*ABS(cn2) - ABS(y-0.5D0*cn2))
  IF (cn1 < qm) THEN
    qm = 1.25D0 * (r - MAX(1.2D0*r, ABS(y-cn2)))
    IF (cn1 < qm) n = 1
  END IF
END IF

qn = pn
a = 4.d-3 * r * r
zz = z
IF (x < 0.0D0) zz = -z

!          CALCULATION OF ZONE OF EXCLUSION OF DEBYE APPROXIMATION

30 nu = CMPLX(qn+fn,cn2)
zn = nu / z
t2 = AIMAG(zn) * AIMAG(zn)
u = 1.0D0 - REAL(zn)
t1 = u * u + t2
u = 1.0D0 + DBLE(zn)
t2 = u * u + t2
u = t1 * t2
v = a * u / (t1*t1 + t2*t2)
IF (u*v*v <= 1.0D0) THEN
  
!          THE ARGUMENT LIES INSIDE THE ZONE OF EXCLUSION
  
  qn = qn + 1.0D0
  IF (n == 0) GO TO 30
  
!          USE MACLAURIN EXPANSION WITH FORWARD RECURRENCE
  
  qn = MIN(pn, REAL(-INT(1.25D0*(r-ABS(cn2)))))
ELSE
  
!          USE BACKWARD RECURRENCE STARTING FROM THE ASYMPTOTIC EXPANSION
  
  qnp1 = qn + 1.0D0
  IF (ABS(qn) < ABS(pn)) THEN
    IF (r >= 17.5D0 + 0.5D0*(qnp1*qnp1 + cn2*cn2)) THEN
      
      nu = CMPLX(qn+fn,cn2)
      CALL cbja(zz,nu,sm1)
      nu = CMPLX(qnp1+fn,cn2)
      CALL cbja(zz,nu,sm2)
      GO TO 40
    END IF
  END IF
  
!          USE BACKWARD RECURRENCE STARTING FROM THE DEBYE APPROXIMATION
  
  nu = CMPLX(qn+fn,cn2)
  CALL cbdb(zz,nu,fn,sm1)
  IF (qn == pn) GO TO 50
  nu = CMPLX(qnp1+fn,cn2)
  CALL cbdb(zz,nu,fn,sm2)
  
  40 nu = CMPLX(qn+fn,cn2)
  tsc = 2.0D0 * nu * sm1 / zz - sm2
  sm2 = sm1
  sm1 = tsc
  qn = qn - 1.0D0
  IF (qn /= pn) GO TO 40
  
  50 w = sm1
  IF (sn < 0.0D0) w = -w
  IF (x >= 0.0D0) RETURN
  
  nu = pi * CMPLX(-cn2,cn1)
  IF (y < 0.0D0) nu = -nu
  w = EXP(nu) * w
  RETURN
END IF

!          USE MACLAURIN EXPANSION WITH FORWARD OR BACKWARD RECURRENCE.

60 m = qn - pn
IF (ABS(m) <= 1) THEN
  nu = CMPLX(cn1,cn2)
  CALL cbjm(z,nu,w)
ELSE
  nu = CMPLX(qn+fn,cn2)
  CALL cbjm(z,nu,w1)
  w0 = 0.25D0 * z * z
  IF (m <= 0) THEN
    
!          FORWARD RECURRENCE
    
    m = ABS(m)
    nu = nu + 1.0D0
    CALL cbjm(z,nu,w)
    DO  i = 2, m
      c = nu * (nu+1.0D0)
      t = (c/w0) * (w-w1)
      w1 = w
      w = t
      nu = nu + 1.0D0
    END DO
  ELSE
    
!          BACKWARD RECURRENCE
    
    nu = nu - 1.0D0
    CALL cbjm(z,nu,w)
    DO  i = 2, m
      c = nu * (nu+1.0D0)
      t = (w0/c) * w1
      w1 = w
      w = w - t
      nu = nu - 1.0D0
    END DO
  END IF
END IF

!          FINAL ASSEMBLY

IF (fn == 0.0D0 .AND. cn2 == 0.0D0) THEN
  k = pn
  IF (k == 0) RETURN
  e = sn / dgamma(real(pn+1.0d0))
  w = e * w * (0.5D0*z) ** k
  RETURN
END IF

s = cnu * LOG(0.5D0*z)
w = EXP(s) * w
IF (rn2 <= 0.81D0) THEN
  w = w * cgam0(cnu)
  RETURN
END IF
CALL dcgama(0,cnu,t)
w = cdiv(w,cnu*t)
RETURN
END SUBROUTINE cbsslj


!-----------------------------------------------------------------------------!



elemental FUNCTION cpabs(x, y) RESULT(fn_val)
!     --------------------------------------
!     EVALUATION OF SQRT(X*X + Y*Y)
!     --------------------------------------
REAL , INTENT(IN)  :: x, y
REAL               :: fn_val

! Local variable
REAL   :: a

IF (ABS(x) > ABS(y)) THEN
  a = y / x
  fn_val = ABS(x) * SQRT(1.0d0 + a*a)
  RETURN
END IF
IF (y /= 0.0d0) THEN
  a = x / y
  fn_val = ABS(y) * SQRT(1.0d0 + a*a)
  RETURN
END IF
fn_val = 0.0d0
RETURN
END FUNCTION cpabs


!-----------------------------------------------------------------------------!



elemental FUNCTION dgamma(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!                EVALUATION OF THE GAMMA FUNCTION FOR
!                     REAL  ARGUMENTS

!                           -----------

!     DGAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL , INTENT(IN) :: a
REAL              :: fn_val

! Local variables
REAL , PARAMETER :: d = 0.41893853320467274178032973640562d0,  &
                        pi = 3.14159265358979323846264338327950d0
REAL  :: s, t, x, w
INTEGER   :: j, n
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
fn_val = 0.0d0
x = a
IF (ABS(a) <= 20.d0) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DGAMMA(A) FOR ABS(A) <= 20
!-----------------------------------------------------------------------
  t = 1.0d0
  n = x
  n = n - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (n < 0) THEN
    GO TO 40
  ELSE IF (n == 0) THEN
    GO TO 30
  END IF

  DO j = 1, n
    x = x - 1.d0
    t = x * t
  END DO
  30 x = x - 1.d0
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0.d0) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1.d0
        t = x * t
      END DO
    END IF
    x = (x+0.5d0) + 0.5d0
    t = x * t
    IF (t == 0.d0) RETURN
  END IF

!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.d-33) THEN
    IF (ABS(t)*HUGE(1.0d0) <= 1.000000001d0) RETURN
    fn_val = 1.d0 / t
    RETURN
  END IF

!     COMPUTE DGAMMA(1 + X) FOR 0 <= X < 1

  60 fn_val = 1.d0 / (1.d0 + dgam1(x))

!     TERMINATION

  IF (a >= 1.d0) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DGAMMA(A) FOR ABS(A) > 20
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.d3) RETURN
IF (a <= 0.d0) THEN
  s = dsin1(a) / pi
  IF (s == 0.d0) RETURN
  x = -a
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

w = dpdel(x)

!     FINAL ASSEMBLY

w = (d+w) + (x-0.5d0) * (LOG(x)-1.d0)
IF (w > dxparg(0)) RETURN
fn_val = EXP(w)
IF (a < 0.d0) fn_val = (1.d0/(fn_val*s)) / x

RETURN
END FUNCTION dgamma

!-----------------------------------------------------------------------------!


elemental FUNCTION cgam0(z) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF 1/GAMMA(1 + Z)  FOR ABS(Z) < 1.0
!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)  :: z
COMPLEX               :: fn_val

COMPLEX   :: w
INTEGER       :: i, k, n
!-----------------------
REAL   :: x, y
REAL, PARAMETER :: a(25) = (/ .577215664901533d0, -.655878071520254d0,  &
        -.420026350340952D-01, .166538611382291d0, -.421977345555443D-01,  &
        -.962197152787697D-02, .721894324666310D-02, -.116516759185907D-02,  &
        -.215241674114951D-03, .128050282388116D-03, -.201348547807882D-04,  &
        -.125049348214267D-0, .113302723198170D-05, -.205633841697761D-0,  &
        .611609510448142D-08, .500200764446922D-08, -.118127457048702D-08,  &
        .104342671169110D-09, .778226343990507D-11, -.369680561864221D-11,  &
        .510037028745448D-12, -.205832605356651D-13, -.534812253942302D-14,  &
        .122677862823826D-14, -.118125930169746D-15 /)
!-----------------------
n = 25
x = REAL(z)
y = AIMAG(z)
IF (x*x+y*y <= 0.04D0) n = 14

k = n
w = a(n)
DO  i = 2, n
  k = k - 1
  w = a(k) + z * w
END DO
fn_val = 1.0D0 + z * w
RETURN
END FUNCTION cgam0


!-----------------------------------------------------------------------------!


elemental FUNCTION cdiv(a,b) RESULT(fn_val)
!-----------------------------------------------------------------------
!              COMPLEX DIVISION A/B WHERE B IS NONZERO
!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)  :: a, b
COMPLEX               :: fn_val


REAL   :: ai, ar, bi, br, d, t
REAL   :: u, v

ar = REAL(a)
ai = AIMAG(a)
br = REAL(b)
bi = AIMAG(b)

IF (ABS(br) >= ABS(bi)) THEN
  t = bi / br
  d = br + t * bi
  u = (ar+ai*t) / d
  v = (ai-ar*t) / d
  fn_val = CMPLX(u,v)
  RETURN
END IF
t = br / bi
d = bi + t * br
u = (ar*t+ai) / d
v = (ai*t-ar) / d
fn_val = CMPLX(u,v)
RETURN
END FUNCTION cdiv



!-----------------------------------------------------------------------------!



elemental FUNCTION dgam1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF 1/GAMMA(1 + X) - 1  FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------

!     THE FOLLOWING ARE THE FIRST 49 COEFFICIENTS OF THE MACLAURIN
!     EXPANSION FOR 1/GAMMA(1 + X) - 1. THE COEFFICIENTS ARE
!     CORRECT TO 40 DIGITS. THE COEFFICIENTS WERE OBTAINED BY
!     ALFRED H. MORRIS JR. (NAVAL SURFACE WARFARE CENTER) AND ARE
!     GIVEN HERE FOR REFERENCE. ONLY THE FIRST 14 COEFFICIENTS ARE
!     USED IN THIS CODE.

!                           -----------

!     DATA A(1)  / .5772156649015328606065120900824024310422D+00/,
!    *     A(2)  /-.6558780715202538810770195151453904812798D+00/,
!    *     A(3)  /-.4200263503409523552900393487542981871139D-01/,
!    *     A(4)  / .1665386113822914895017007951021052357178D+00/,
!    *     A(5)  /-.4219773455554433674820830128918739130165D-01/,
!    *     A(6)  /-.9621971527876973562114921672348198975363D-02/,
!    *     A(7)  / .7218943246663099542395010340446572709905D-02/,
!    *     A(8)  /-.1165167591859065112113971084018388666809D-02/,
!    *     A(9)  /-.2152416741149509728157299630536478064782D-03/,
!    *     A(10) / .1280502823881161861531986263281643233949D-03/
!     DATA A(11) /-.2013485478078823865568939142102181838229D-04/,
!    *     A(12) /-.1250493482142670657345359473833092242323D-05/,
!    *     A(13) / .1133027231981695882374129620330744943324D-05/,
!    *     A(14) /-.2056338416977607103450154130020572836513D-06/,
!    *     A(15) / .6116095104481415817862498682855342867276D-08/,
!    *     A(16) / .5002007644469222930055665048059991303045D-08/,
!    *     A(17) /-.1181274570487020144588126565436505577739D-08/,
!    *     A(18) / .1043426711691100510491540332312250191401D-09/,
!    *     A(19) / .7782263439905071254049937311360777226068D-11/,
!    *     A(20) /-.3696805618642205708187815878085766236571D-11/
!     DATA A(21) / .5100370287454475979015481322863231802727D-12/,
!    *     A(22) /-.2058326053566506783222429544855237419746D-13/,
!    *     A(23) /-.5348122539423017982370017318727939948990D-14/,
!    *     A(24) / .1226778628238260790158893846622422428165D-14/,
!    *     A(25) /-.1181259301697458769513764586842297831212D-15/,
!    *     A(26) / .1186692254751600332579777242928674071088D-17/,
!    *     A(27) / .1412380655318031781555803947566709037086D-17/,
!    *     A(28) /-.2298745684435370206592478580633699260285D-18/,
!    *     A(29) / .1714406321927337433383963370267257066813D-19/,
!    *     A(30) / .1337351730493693114864781395122268022875D-21/
!     DATA A(31) /-.2054233551766672789325025351355733796682D-21/,
!    *     A(32) / .2736030048607999844831509904330982014865D-22/,
!    *     A(33) /-.1732356445910516639057428451564779799070D-23/,
!    *     A(34) /-.2360619024499287287343450735427531007926D-25/,
!    *     A(35) / .1864982941717294430718413161878666898946D-25/,
!    *     A(36) /-.2218095624207197204399716913626860379732D-26/,
!    *     A(37) / .1297781974947993668824414486330594165619D-27/,
!    *     A(38) / .1180697474966528406222745415509971518560D-29/,
!    *     A(39) /-.1124584349277088090293654674261439512119D-29/,
!    *     A(40) / .1277085175140866203990206677751124647749D-30/
!     DATA A(41) /-.7391451169615140823461289330108552823711D-32/,
!    *     A(42) / .1134750257554215760954165259469306393009D-34/,
!    *     A(43) / .4639134641058722029944804907952228463058D-34/,
!    *     A(44) /-.5347336818439198875077418196709893320905D-35/,
!    *     A(45) / .3207995923613352622861237279082794391090D-36/,
!    *     A(46) /-.4445829736550756882101590352124643637401D-38/,
!    *     A(47) /-.1311174518881988712901058494389922190237D-38/,
!    *     A(48) / .1647033352543813886818259327906394145400D-39/,
!    *     A(49) /-.1056233178503581218600561071538285049997D-40/

!                           -----------

!     C = A(1) - 1 IS ALSO FREQUENTLY NEEDED. C HAS THE VALUE ...

!     DATA C /-.4227843350984671393934879099175975689578D+00/

!-----------------------------------------------------------------------
REAL , INTENT(IN) :: x
REAL              :: fn_val

! Local variables
REAL  :: d, t, w, z
REAL , PARAMETER :: a0 = .611609510448141581788D-08, a1  &
        = .624730830116465516210D-08, b1 = .203610414066806987300D+00, b2  &
        = .266205348428949217746D-01, b3 = .493944979382446875238D-03, b4  &
        = -.851419432440314906588D-05, b5 = -.643045481779353022248D-05, b6  &
        = .992641840672773722196D-06, b7 = -.607761895722825260739D-07, b8  &
        = .195755836614639731882D-09
REAL , PARAMETER :: p0 = .6116095104481415817861D-08, p1  &
        = .6871674113067198736152D-08, p2 = .6820161668496170657, p3  &
        = .4686843322948848031080D-10, p4 = .1572833027710446286995D-11, p5  &
        = -.1249441572276366213222D-12, p6 = .4343529937408594255178D-14, q1  &
        = .3056961078365221025009D+00, q2 = .5464213086042296536016D-01, q3  &
        = .4956830093825887312, q4 = .2692369466186361192876D-03
REAL , PARAMETER :: c = -.422784335098467139393487909917598D+00, c0  &
        = .577215664901532860606512090082402D+00, c1  &
        = -.655878071520253881077019515145390D+00, c2  &
        = -.420026350340952355290039348754298D-01, c3  &
        = .166538611382291489501700795102105D+00, c4  &
        = -.421977345555443367482083012891874D-01, c5  &
        = -.962197152787697356211492167234820D-02, c6  &
        = .721894324666309954239501034044657D-02, c7  &
        = -.116516759185906511211397108401839D-02, c8  &
        = -.215241674114950972815729963053648D-03, c9  &
        = .128050282388116186153198626328164D-03, c10  &
        = -.201348547807882386556893914210218D-04, c11  &
        = -.125049348214267065734535947383309D-05, c12  &
        = .113302723198169588237412962033074D-05, c13  &
        = -.205633841697760710345015413002057D-06
!----------------------------
t = x
d = x - 0.5d0
IF (d > 0.d0) t = d - 0.5d0
IF (t < 0.0d0) THEN
  GO TO 30
ELSE IF (t > 0.0d0) THEN
  GO TO 20
END IF

fn_val = 0.d0
RETURN
!------------

!                CASE WHEN 0 < T <= 0.5

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
20 w = ((((((p6*t + p5)*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) /   &
       ((((q4*t+q3)*t + q2)*t + q1)*t + 1.d0)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c0

IF (d <= 0.d0) THEN
  fn_val = x * z
  RETURN
END IF
fn_val = (t/x) * ((z-0.5d0)-0.5d0)
RETURN
!------------

!                CASE WHEN -0.5 <= T < 0

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
30 w = (a1*t + a0) / ((((((((b8*t + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t+1.d0)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c

IF (d <= 0.d0) THEN
  fn_val = x * ((z+0.5d0)+0.5d0)
  RETURN
END IF
fn_val = t * z / x
RETURN
END FUNCTION dgam1


!-----------------------------------------------------------------------------!


elemental FUNCTION dsin1(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!                REAL  EVALUATION OF SIN(PI*X)

!                             --------------

!     THE EXPANSION FOR SIN(PI*A) (ABS(A) <= PI/4) USING A1,...,A13
!     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!     THE EXPANSION FOR COS(PI*A) (ABS(A) <= PI/4) USING B1,...,B13
!     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!-----------------------------------------------------------------------
REAL , INTENT(IN)  :: x
REAL               :: fn_val

! Local variables
REAL             :: a, t, w
REAL , PARAMETER :: pi = 3.141592653589793238462643383279502884197D+00
REAL , PARAMETER :: a1 = -.1028083791780141522795259479153765743002D+00,   &
      a2  = .3170868848763100170457042079710451905600D-02,   &
      a3  = -.4657026956105571623449026167864697920000D-04,  &
      a4  = .3989844942879455643410226655783424000000D-06,   &
      a5  = -.2237397227721999776371894030796800000000D-08,  &
      a6  = .8847045483056962709715066675200000000000D-11,   &
      a7  = -.2598715447506450292885585920000000000000D-13,  &
      a8  = .5893449774331011070033920000000000000000D-16 ,  &
      a9  = -.1062975472045522550784000000000000000000D-18,   &
      a10 = .1561182648301780992000000000000000000000D-21,    &
      a11 = -.1903193516670976000000000000000000000000D-24,   &
      a12 = .1956617650176000000000000000000000000000D-27,    &
      a13 = -.1711276032000000000000000000000000000000D-30
REAL , PARAMETER :: b1 = -.3084251375340424568385778437461297229882D+00, &
      b2  = .1585434424381550085228521039855226435920D-01,   &
      b3  = -.3259918869273900136414318317506279360000D-03,  &
      b4  = .3590860448591510079069203991239232000000D-05,   &
      b5  = -.2461136950494199754009084061808640000000D-07,  &
      b6  = .1150115912797405152263195572224000000000D-09,   &
      b7  = -.3898073171259675439899172864000000000000D-12,  &
      b8  = .1001886461636271969091584000000000000000D-14,   &
      b9  = -.2019653396886572027084800000000000000000D-17,  &
      b10 = .3278483561466560512000000000000000000000D-20,   &
      b11 = -.4377345082051788800000000000000000000000D-23,  &
      b12 = .4891532381388800000000000000000000000000D-26,   &
      b13 = -.4617089843200000000000000000000000000000D-29
INTEGER  :: max, n
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
t = MAX
IF (a >= t) THEN
  fn_val = 0.0d0
  RETURN
END IF

n = a
t = n
a = a - t
IF (a <= 0.75d0) THEN
  IF (a < 0.25d0) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25d0 + (0.25d0-a)
  t = 16.d0 * a * a
  fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t  &
           + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t +  &
           0.5d0) + 0.5d0
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25d0 + (0.75d0-a)
10 t = 16.d0 * a * a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t + a8)*t + a7)*t  &
    + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5d0) + 0.5d0
fn_val = pi * a * w

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dsin1

!-----------------------------------------------------------------------------!


elemental FUNCTION dpdel(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF THE FUNCTION DEL(X) FOR  X >= 10  WHERE
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X)

!                         --------

!     THE SERIES FOR DPDEL ON THE INTERVAL 0.0 TO 1.0 DERIVED BY
!     A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY
!     OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!-----------------------------------------------------------------------
REAL , INTENT(IN) :: x
REAL              :: fn_val

! Local variables
REAL , PARAMETER :: a(15) = (/ .833333333333333333333333333333D-01,  &
        -.277777777777777777777777752282D-04,  &
        .793650793650793650791732130419D-07,  &
        -.595238095238095232389839236182D-09,  &
        .841750841750832853294451671990D-11,  &
        -.191752691751854612334149171243D-12,  &
        .641025640510325475730918472625D-14,  &
        -.295506514125338232839867823991D-15,  &
        .179643716359402238723287696452D-16,  &
        -.139228964661627791231203060395D-17,  &
        .133802855014020915603275339093D-18,  &
        -.154246009867966094273710216533D-19,  &
        .197701992980957427278370133333D-20,  &
        -.234065664793997056856992426667D-21,  &
        .171348014966398575409015466667D-22 /)
REAL  :: t, w
INTEGER   :: i, k
!-----------------------------------------------------------------------
t = (10.d0/x) ** 2
w = a(15)
DO i = 1, 14
  k = 15 - i
  w = t * w + a(k)
END DO
fn_val = w / x
RETURN
END FUNCTION dpdel


elemental FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL            :: fn_val

! Local variable
REAL     :: one     !changed

one=1.0d0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg

!-----------------------------------------------------------------------------!



elemental SUBROUTINE cbja(cz,cnu,w)
!-----------------------------------------------------------------------
!        COMPUTATION OF J(NU,Z) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)   :: cz
COMPLEX , INTENT(IN)   :: cnu
COMPLEX , INTENT(OUT)  :: w

! Local variables
REAL      :: eps, inu, m
COMPLEX   :: a, a1, arg, e, eta, nu, p, q, t, z, zr, zz
!--------------------------
REAL  :: r, rnu, tol, u, v
REAL  :: x, y
INTEGER   :: i, ind

!--------------------------
!     PIHALF = PI/2
!     C = 2*PI**(-1/2)
!--------------------------
REAL , PARAMETER    :: pihalf = 1.5707963267949d0, c = 1.12837916709551d0
COMPLEX , PARAMETER :: j = (0.0d0, 1.0d0)
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0 .

eps = EPSILON(0.0d0)

!--------------------------
z = cz
x = REAL(z)
y = AIMAG(z)
nu = cnu
ind = 0
IF (ABS(x) <= 1.d-2*ABS(y)) THEN
  IF (AIMAG(nu) < 0.0D0 .AND. ABS(REAL(nu)) < 1.d-2*ABS(AIMAG(nu))) THEN
    ind = 1
    nu = CONJG(nu)
    z = CONJG(z)
    y = -y
  END IF
END IF

IF (x < -1.d-2*y) z = -z
zz = z + z
CALL dcrec(REAL(zz),AIMAG(zz),u,v)
zr = CMPLX(u,v)
eta = -zr * zr

p = (0.0D0,0.0D0)
q = (0.0D0,0.0D0)
a1 = nu * nu - 0.25D0
a = a1
t = a1
m = 1.0D0
tol = eps * anorm(a1)
DO  i = 1, 16
  a = a - 2.0D0 * m
  m = m + 1.0D0
  t = t * a * eta / m
  p = p + t
  a = a - 2.0D0 * m
  m = m + 1.0D0
  t = t * a / m
  q = q + t
  IF (anorm(t) <= tol) GO TO 20
END DO

20 p = p + 1.0D0
q = (q+a1) * zr
w = z - pihalf * nu
IF (ABS(AIMAG(w)) <= 1.0D0) THEN
  arg = w - 0.5D0 * pihalf
  w = c * SQRT(zr) * (p*COS(arg) - q*SIN(arg))
ELSE
  e = EXP(-j*w)
  t = q - j * p
  IF (AIMAG(z) > 0.0D0 .AND. REAL(z) <= 1.d-2*AIMAG(z).AND.  &
      ABS(REAL(nu)) < 1.d-2*AIMAG(nu)) t = 0.5D0 * t
  CALL dcrec(REAL(e),AIMAG(e),u,v)
  w = 0.5D0 * c * SQRT(j*zr) * ((p-j*q)*e + t*CMPLX(u,v))
END IF

IF (x < -1.d-2*y) THEN
  IF (y < 0.0D0) nu = -nu
  
!     COMPUTATION OF EXP(I*PI*NU)
  
  rnu = REAL(nu)
  inu = AIMAG(nu)
  r = EXP(-2.0D0*pihalf*inu)
  u = r * dcos1(rnu)
  v = r * dsin1(rnu)
  w = w * CMPLX(u,v)
END IF

IF (ind /= 0) w = CONJG(w)
RETURN
END SUBROUTINE cbja

!-----------------------------------------------------------------------------!


elemental FUNCTION anorm(z) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

COMPLEX , INTENT(IN)  :: z
REAL                  :: fn_val

fn_val = MAX( ABS( REAL(z)), ABS(AIMAG(z) ) )
RETURN
END FUNCTION anorm

!-----------------------------------------------------------------------------!


elemental FUNCTION dcos1 (x) RESULT(fn_val)
 
!-----------------------------------------------------------------------

!                REAL  EVALUATION OF COS(PI*X)

!                             --------------

!     THE EXPANSION FOR SIN(PI*A) (ABS(A) .LE. PI/4) USING A1,...,A13
!     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!     THE EXPANSION FOR COS(PI*A) (ABS(A) .LE. PI/4) USING B1,...,B13
!     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!-----------------------------------------------------------------------

REAL , INTENT(IN)  :: x
REAL               :: fn_val

REAL   :: a, t, w
INTEGER    :: MAX, n
!------------------------
REAL , PARAMETER  :: pi = 3.141592653589793238462643383279502884197d0
!------------------------
REAL , PARAMETER  :: &
    a1  = -.1028083791780141522795259479153765743002D+00,  &
    a2  =  .3170868848763100170457042079710451905600D-02,  &
    a3  = -.4657026956105571623449026167864697920000D-04,  &
    a4  =  .3989844942879455643410226655783424000000D-06,  &
    a5  = -.2237397227721999776371894030796800000000D-08,  &
    a6  =  .8847045483056962709715066675200000000000D-11,  &
    a7  = -.2598715447506450292885585920000000000000D-13,  &
    a8  =  .5893449774331011070033920000000000000000D-16,  &
    a9  = -.1062975472045522550784000000000000000000D-18,  &
    a10 =  .1561182648301780992000000000000000000000D-21,  &
    a11 = -.1903193516670976000000000000000000000000D-24,  &
    a12 =  .1956617650176000000000000000000000000000D-27,  &
    a13 = -.1711276032000000000000000000000000000000D-30
!------------------------
REAL , PARAMETER  :: &
    b1  = -.3084251375340424568385778437461297229882D+00,  &
    b2  =  .1585434424381550085228521039855226435920D-01,  &
    b3  = -.3259918869273900136414318317506279360000D-03,  &
    b4  =  .3590860448591510079069203991239232000000D-05,  &
    b5  = -.2461136950494199754009084061808640000000D-07,  &
    b6  =  .1150115912797405152263195572224000000000D-09,  &
    b7  = -.3898073171259675439899172864000000000000D-12,  &
    b8  =  .1001886461636271969091584000000000000000D-14,  &
    b9  = -.2019653396886572027084800000000000000000D-17,  &
    b10 =  .3278483561466560512000000000000000000000D-20,  &
    b11 = -.4377345082051788800000000000000000000000D-23,  &
    b12 =  .4891532381388800000000000000000000000000D-26,  &
    b13 = -.4617089843200000000000000000000000000000D-29
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

MAX = HUGE(0)

!------------------------
a = ABS(x)
t = MAX
IF (a < t) GO TO 10
fn_val = 1.d0
RETURN

10 n = a
t = n
a = a - t
IF (a > 0.75D0) GO TO 20
IF (a < 0.25D0) GO TO 21

!                    0.25 .LE. A .LE. 0.75

a = 0.25D0 + (0.25D0 - a)
t = 16.d0*a*a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t +  &
    a8)*t + a7)*t + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5D0) + 0.5D0
fn_val = pi*a*w
GO TO 30

!                 A .LT. 0.25  OR  A .GT. 0.75

20 a = 0.25D0 + (0.75D0 - a)
n = n - 1
21 t = 16.d0*a*a
fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t + &
         b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t + 0.5D0) + 0.5D0

!                        TERMINATION

30 IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dcos1

!-----------------------------------------------------------------------------!


elemental SUBROUTINE cbjm(z,cnu,w)
!-----------------------------------------------------------------------

!       COMPUTATION OF  (Z/2)**(-CNU) * GAMMA(CNU + 1) * J(CNU,Z)

!                           -----------------

!     THE MACLAURIN EXPANSION IS USED.  IT IS ASSUMED THAT CNU IS NOT
!     A NEGATIVE INTEGER.

!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)   :: z
COMPLEX , INTENT(IN)   :: cnu
COMPLEX , INTENT(OUT)  :: w

COMPLEX   :: nu, nup1, p, s, sn, t, ti
!--------------------------
REAL   :: a, a0, eps, inu, m, rnu
INTEGER    :: i, imin, k, km1, km2

!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0 .

eps = EPSILON(0.0d0)

!--------------------------
s = -0.25D0 * (z*z)
nu = cnu
rnu = REAL(nu)
inu = AIMAG(nu)
a = 0.5D0 + (0.5D0+rnu)
nup1 = CMPLX(a,inu)

IF (a > 0.0D0) THEN
  m = 1.0D0
  t = s / nup1
  w = 1.0D0 + t
ELSE
  
!     ADD 1.0 AND THE FIRST K-1 TERMS
  
  k = INT(-a) + 2
  km1 = k - 1
  w = (1.0D0,0.0D0)
  t = w
  DO  i = 1, km1
    m = i
    t = t * (s/(m*(nu+m)))
    w = w + t
    IF (anorm(t) <= eps*anorm(w)) GO TO 20
  END DO
  GO TO 50
  
!     CHECK IF THE (K-1)-ST AND K-TH TERMS CAN BE IGNORED.
!     IF SO THEN THE SUMMATION IS COMPLETE.
  
  20 IF (i /= km1) THEN
    imin = i + 1
    IF (imin < k-5) THEN
      ti = t
      
      m = km1
      t = s / (nu+m)
      a0 = anorm(t) / m
      t = t * (s/(nu+(m+1.0D0)))
      a = anorm(t) / (m*(m+1.0D0))
      a = MAX(a,a0)
      
      t = (1.0D0,0.0D0)
      km2 = k - 2
      DO  i = imin, km2
        m = i
        t = t * (s/(m*(nu+m)))
        IF (a*anorm(t) < 0.5D0) RETURN
      END DO
      t = t * ti
      imin = km2
    END IF
    
!     ADD THE (K-1)-ST TERM
    
    a = 1.0D0
    p = (1.0D0,0.0D0)
    sn = p
    DO  i = imin, km1
      m = i
      a = a * m
      p = p * (nu+m)
      sn = s * sn
    END DO
    t = t * (cdiv(sn,p)/a)
    w = w + t
  END IF
END IF

!     ADD THE REMAINING TERMS

50 m = m + 1.0D0
t = t * (s/(m*(nu+m)))
w = w + t
IF (anorm(t) > eps*anorm(w)) GO TO 50

RETURN
END SUBROUTINE cbjm

!-----------------------------------------------------------------------------!


elemental SUBROUTINE dcgama(mo, z, w)
!-----------------------------------------------------------------------

!        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

!                        ---------------

!     MO IS AN INTEGER.  Z AND W ARE INTERPRETED AS REAL 
!     COMPLEX NUMBERS.  IT IS ASSUMED THAT Z(1) AND Z(2) ARE THE REAL
!     AND IMAGINARY PARTS OF THE COMPLEX NUMBER Z, AND THAT W(1) AND
!     W(2) ARE THE REAL AND IMAGINARY PARTS OF W.

!                 W = GAMMA(Z)       IF MO = 0
!                 W = LN(GAMMA(Z))   OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)        :: mo
COMPLEX , INTENT(IN)   :: z
COMPLEX , INTENT(OUT)  :: w

! Local variables
REAL , PARAMETER :: c0(30)  &
        = (/ .8333333333333333333333333333333333333333D-01,  &
        -.2777777777777777777777777777777777777778D-02,  &
         .7936507936507936507936507936507936507937D-03,  &
        -.5952380952380952380952380952380952380952D-03,  &
         .8417508417508417508417508417508417508418D-03,  &
        -.1917526917526917526917526917526917526918D-02,  &
         .6410256410256410256410256410256410256410D-02,  &
        -.2955065359477124183006535947712418300654D-01,  &
         .1796443723688305731649384900158893966944D+00,  &
        -.1392432216905901116427432216905901116427D+01,  &
         .1340286404416839199447895100069013112491D+02,  &
        -.1568482846260020173063651324520889738281D+03,  &
         .2193103333333333333333333333333333333333D+04,  &
        -.3610877125372498935717326521924223073648D+05,  &
         .6914722688513130671083952507756734675533D+06,  &
        -.1523822153940741619228336495888678051866D+08,  &
         .3829007513914141414141414141414141414141D+09,  &
        -.1088226603578439108901514916552510537473D+11,  &
         .3473202837650022522522522522522522522523D+12,  &
        -.1236960214226927445425171034927132488108D+14,  &
         .4887880647930793350758151625180229021085D+15,  &
        -.2132033396091937389697505898213683855747D+17,  &
         .1021775296525700077565287628053585500394D+19,  &
        -.5357547217330020361082770919196920448485D+20,  &
         .3061578263704883415043151051329622758194D+22,  &
        -.1899991742639920405029371429306942902947D+24,  &
         .1276337403382883414923495137769782597654D+26,  &
        -.9252847176120416307230242348347622779519D+27,  &
         .7218822595185610297836050187301637922490D+29,  &
        -.6045183405995856967743148238754547286066D+31 /),  &
        dlpi = 1.144729885849400174143427351353058711647d0,  &
        hl2p =  .9189385332046727417803297364056176398614d0,  &
        pi = 3.141592653589793238462643383279502884197d0,  &
        pi2 = 6.283185307179586476925286766559005768394d0
REAL  :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, q1, q2, s, sn,  &
             s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
INTEGER   :: j, k, max, n, nm1
!---------------------------
!     DLPI = LOG(PI)
!     HL2P = 0.5 * LOG(2*PI)
!---------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST NUMBER SUCH
!            THAT  1.d0 + EPS > 1.d0.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0d0)

!---------------------------
x = REAL(z)
y = AIMAG(z)
IF (x < 0.d0) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi * y
  et = EXP(t)
  e2t = et * et

!     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

  a1 = 0.5d0 * (1.d0+e2t)
  t2 = t + t
  IF (t2 >= -0.15d0) THEN
    a2 = -0.5d0 * drexp(t2)
  ELSE
    a2 = 0.5d0 * (0.5d0+(0.5d0-e2t))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X)

  u = MAX
  IF (ABS(x) >= MIN(u,1.d0/eps)) GO TO 80
  k = ABS(x)
  u = x + k
  k = MOD(k,2)
  IF (u <= -0.5d0) THEN
    u = 0.5d0 + (0.5d0+u)
    k = k + 1
  END IF
  u = pi * u
  sn = SIN(u)
  cn = COS(u)
  IF (k == 1) THEN
    sn = -sn
    cn = -cn
  END IF

!     SET  H1 + H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))

  a1 = sn * a1
  a2 = cn * a2
  a = a1 * a1 + a2 * a2
  IF (a == 0.d0) GO TO 80
  IF (mo == 0) THEN

    h1 = a1 / a
    h2 = -a2 / a
    c = pi * et
    h1 = c * h1
    h2 = c * h2
  ELSE

    h1 = (dlpi+t) - 0.5d0 * LOG(a)
    h2 = -ATAN2(a2,a1)
  END IF
  IF (AIMAG(z) >= 0.d0) THEN
    x = 1.0 - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0 - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
w1 = 0.d0
w2 = 0.d0
n = 0
t = x
y2 = y * y
a = t * t + y2
cut = 225.d0
IF (eps > 1.d-30) cut = 144.d0
IF (eps > 1.d-20) cut = 64.d0
IF (a < cut) THEN
  IF (a == 0.d0) GO TO 80
  10 n = n + 1
  t = t + 1.d0
  a = t * t + y2
  IF (a < cut) GO TO 10

!     LET S1 + S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)

  u1 = (x*t+y2) / a
  u2 = y / a
  s1 = u1
  s2 = n * u2
  IF (n >= 2) THEN
    u = t / a
    nm1 = n - 1
    DO j = 1, nm1
      v1 = u1 + j * u
      v2 = (n-j) * u2
      c = s1 * v1 - s2 * v2
      d = s1 * v2 + s2 * v1
      s1 = c
      s2 = d
    END DO
  END IF

!     SET  W1 + W2*I = LOG(S1 + S2*I)  WHEN MO IS NONZERO

  s = s1 * s1 + s2 * s2
  IF (mo /= 0) THEN
    w1 = 0.5d0 * LOG(s)
    w2 = ATAN2(s2,s1)
  END IF
END IF

!     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

t1 = 0.5d0 * LOG(a) - 1.d0
t2 = ATAN2(y,t)
u = x - 0.5d0
v1 = (u*t1-0.5d0) - y * t2
v2 = u * t2 + y * t1

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

u1 = t / a
u2 = -y / a
q1 = u1 * u1 - u2 * u2
q2 = 2.d0 * u1 * u2
a1 = 0.d0
a2 = 0.d0
DO j = 1, 30
  t1 = a1
  t2 = a2
  a1 = a1 + c0(j) * u1
  a2 = a2 + c0(j) * u2
  IF (a1 == t1) THEN
    IF (a2 == t2) GO TO 40
  END IF
  t1 = u1 * q1 - u2 * q2
  t2 = u1 * q2 + u2 * q1
  u1 = t1
  u2 = t2
END DO
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
40 w1 = (((a1+hl2p)-w1)+v1) - n
w2 = (a2-w2) + v2
IF (REAL(z) < 0.0d0) GO TO 60
IF (mo == 0) THEN

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

  a = EXP(w1)
  w1 = a * COS(w2)
  w2 = a * SIN(w2)
  IF (n == 0) GO TO 70
  c = (s1*w1+s2*w2) / s
  d = (s1*w2-s2*w1) / s
  w1 = c
  w2 = d
  GO TO 70
END IF

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
!     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

50 IF (w2 <= pi) THEN
  k = 0.5d0 - w2 / pi2
  w2 = w2 + pi2 * k
  GO TO 70
END IF
k = w2 / pi2 - 0.5d0
u = k + 1
w2 = w2 - pi2 * u
IF (w2 <= -pi) w2 = pi
GO TO 70

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

60 IF (mo /= 0) THEN
  w1 = h1 - w1
  w2 = h2 - w2
  GO TO 50
END IF

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO = 0

a = EXP(-w1)
t1 = a * COS(-w2)
t2 = a * SIN(-w2)
w1 = h1 * t1 - h2 * t2
w2 = h1 * t2 + h2 * t1
IF (n /= 0) THEN
  c = w1 * s1 - w2 * s2
  d = w1 * s2 + w2 * s1
  w1 = c
  w2 = d
END IF

!     TERMINATION

70 w = CMPLX(w1, w2)
RETURN
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
80 w = CMPLX(0.0d0, 0.0d0)
RETURN
END SUBROUTINE dcgama

!-----------------------------------------------------------------------------!


elemental FUNCTION drexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL , INTENT(IN)  :: x
REAL               :: fn_val

! Local variables
REAL  :: e, w, z
REAL, parameter  :: a0 = .248015873015873015873016D-04,   &
    a1 = -.344452080605731005808147D-05, a2 = .206664230430046597475413D-06,  &
    a3 = -.447300111094328162971036D-08, a4 = .114734027080634968083920D-11,  &
    b1 = -.249994190011341852652396D+00, b2 = .249987228833107957725728D-01,  &
    b3 = -.119037506846942249362528D-02, b4 = .228908693387350391768682D-04
REAL, parameter :: c1 = .1666666666666666666666666666666667D+00,   &
             c2 = .4166666666666666666666666666666667D-01,   &
             c3 = .8333333333333333333333333333333333D-02,   &
             c4 = .1388888888888888888888888888888889D-02,   &
             c5 = .1984126984126984126984126984126984D-03
!---------------------------
IF (ABS(x) <= 0.15d0) THEN

!        Z IS A MINIMAX APPROXIMATION OF THE SERIES

!                C6 + C7*X + C8*X**2 + ....

!        THIS APPROXIMATION IS ACCURATE TO WITHIN
!        1 UNIT OF THE 23-RD SIGNIFICANT DIGIT.
!        THE RESULTING VALUE FOR W IS ACCURATE TO
!        WITHIN 1 UNIT OF THE 33-RD SIGNIFICANT DIGIT.

  z = ((((a4*x + a3)*x + a2)*x + a1)*x + a0) /  &
      ((((b4*x + b3)*x + b2)*x + b1)*x + 1.d0)
  w = ((((((z*x + c5)*x + c4)*x + c3)*x + c2)*x + c1)*x + 0.5d0)*x + 1.d0
  fn_val = x * w
  RETURN
END IF

IF (x >= 0.0d0) THEN
  e = EXP(x)
  fn_val = e * (0.5d0 + (0.5d0 - 1.0d0/e))
  RETURN
END IF
IF (x >= -77.d0) THEN
  fn_val = (EXP(x)-0.5d0) - 0.5d0
  RETURN
END IF
fn_val = -1.d0
RETURN
END FUNCTION drexp

!-----------------------------------------------------------------------------!


elemental SUBROUTINE cbdb(cz,cnu,fn,w)
!-----------------------------------------------------------------------

!         CALCULATION OF J   (CZ) BY THE DEBYE APPROXIMATION
!                         CNU
!                         ------------------

!     IT IS ASSUMED THAT REAL(CZ) .GE. 0 AND THAT REAL(CNU) = FN + K
!     WHERE K IS AN INTEGER.

!-----------------------------------------------------------------------

COMPLEX , INTENT(IN)   :: cz, cnu
REAL , INTENT(IN)      :: fn
COMPLEX , INTENT(OUT)  :: w

! Local variables
REAL      :: is, inu, izn
COMPLEX   :: c1, c2, eta, nu, p, p1, q, r, s, s1, s2, sm, t, z, zn
!----------------------
!     C = 1/SQRT(2*PI)
!     BND = PI/3
!----------------------
REAL , PARAMETER  :: c = .398942280401433d0, pi = 3.14159265358979d0,  &
                         pi2 = 6.28318530717959d0, bnd = 1.04719755119660d0
COMPLEX , PARAMETER :: j = (0.0, 1.0)
REAL   :: alpha, am, aq, ar
REAL   :: phi, sgn, theta
REAL   :: u, v, x, y
INTEGER    :: ind, k, l, m

!----------------------
!             COEFFICIENTS OF THE FIRST 16 POLYNOMIALS
!                   IN THE DEBYE APPROXIMATION
!----------------------

REAL, parameter   :: a(136) = (/ 1.0d0, -.208333333333333d0, .125000000000000d0, .334201388888889d0, &
  -.401041666666667d0, .703125000000000D-01,-.102581259645062D+01, .184646267361111D+01, &
  -.891210937500000d0, .732421875000000D-01, .466958442342625D+01,-.112070026162230D+02, &
   .878912353515625D+01,-.236408691406250D+01, .112152099609375d0,-.282120725582002D+02, &
   .846362176746007D+02,-.918182415432400D+02, .425349987453885D+02,-.736879435947963D+01, &
   .227108001708984d0, .212570130039217D+03,-.765252468141182D+03, .105999045252800D+04, &
  -.699579627376133D+03, .218190511744212D+03,-.264914304869516D+02, .572501420974731d0, &
  -.191945766231841D+04, .806172218173731D+04,-.135865500064341D+05, .116553933368645D+05, &
  -.530564697861340D+04, .120090291321635D+04,-.108090919788395D+03, .172772750258446D+01, &
   .202042913309661D+05,-.969805983886375D+05, .192547001232532D+06,-.203400177280416D+06, &
   .122200464983017D+06,-.411926549688976D+05, .710951430248936D+04,-.493915304773088D+03, &
   .607404200127348D+01,-.242919187900551D+06, .131176361466298D+07,-.299801591853811D+07, &
   .376327129765640D+07,-.281356322658653D+07, .126836527332162D+07,-.331645172484564D+06, &
   .452187689813627D+05,-.249983048181121D+04, .243805296995561D+02, .328446985307204D+07, &
  -.197068191184322D+08, .509526024926646D+08,-.741051482115327D+08, .663445122747290D+08, &
  -.375671766607634D+08, .132887671664218D+08,-.278561812808645D+07, .308186404612662D+06, &
  -.138860897537170D+05, .110017140269247D+03,-.493292536645100D+08, .325573074185766D+09, &
  -.939462359681578D+09, .155359689957058D+10,-.162108055210834D+10, .110684281682301D+10, &
  -.495889784275030D+09, .142062907797533D+09,-.244740627257387D+08, .224376817792245D+07, &
  -.840054336030241D+05, .551335896122021D+03, .814789096118312D+09,-.586648149205185D+10, &
   .186882075092958D+11,-.346320433881588D+11, .412801855797540D+11,-.330265997498007D+11, &
   .179542137311556D+11,-.656329379261928D+10, .155927986487926D+10,-.225105661889415D+09, &
   .173951075539782D+08,-.549842327572289D+06, .303809051092238D+04,-.146792612476956D+11, &
   .114498237732026D+12,-.399096175224466D+12, .819218669548577D+12,-.109837515608122D+13, &
   .100815810686538D+13,-.645364869245377D+12, .287900649906151D+12,-.878670721780233D+11, &
   .176347306068350D+11,-.216716498322380D+10, .143157876718889D+09,-.387183344257261D+07, &
   .182577554742932D+05, .286464035717679D+12,-.240629790002850D+13, .910934118523990D+13, &
  -.205168994109344D+14, .305651255199353D+14,-.316670885847852D+14, .233483640445818D+14, &
  -.123204913055983D+14, .461272578084913D+13,-.119655288019618D+13, .205914503232410D+12, &
  -.218229277575292D+11, .124700929351271D+10,-.291883881222208D+08, .118838426256783D+06, &
  -.601972341723401D+13, .541775107551060D+14,-.221349638702525D+15, .542739664987660D+15, &
  -.889496939881026D+15, .102695519608276D+16,-.857461032982895D+15, .523054882578445D+15, &
  -.232604831188940D+15, .743731229086791D+14,-.166348247248925D+14, .248500092803409D+13, &
  -.229619372968246D+12, .114657548994482D+11,-.234557963522252D+09, .832859304016289D+06 /)

z = cz
nu = cnu
inu = AIMAG(cnu)
IF (inu < 0.0D0) THEN
  z = CONJG(z)
  nu = CONJG(nu)
END IF
x = REAL(z)
y = AIMAG(z)

!          TANH(GAMMA) = SQRT(1 - (Z/NU)**2) = W/NU
!          T = EXP(NU*(TANH(GAMMA) - GAMMA))

zn = z / nu
izn = AIMAG(zn)
IF (ABS(izn) <= 0.1D0*ABS(REAL(zn))) THEN
  
  s = (1.0D0-zn) * (1.0D0+zn)
  eta = 1.0D0 / s
  q = SQRT(s)
  s = 1.0D0 / (nu*q)
  t = zn / (1.0D0 + q)
  t = EXP(nu*(q + LOG(t)))
ELSE
  
  s = (nu-z) * (nu+z)
  eta = (nu*nu) / s
  w = SQRT(s)
  q = w / nu
  IF (REAL(q) < 0.0D0) w = -w
  s = 1.0D0 / w
  t = z / (nu+w)
  t = EXP(w + nu*LOG(t))
END IF

is = AIMAG(s)
r = SQRT(s)
c1 = r * t
ar = REAL(r) * REAL(r) + AIMAG(r) * AIMAG(r)
aq = -1.0D0 / (REAL(q)*REAL(q) + AIMAG(q)*AIMAG(q))

phi = ATAN2(y,x) / 3.0D0
q = nu - z
theta = ATAN2(AIMAG(q),REAL(q)) - phi
ind = 0
IF (ABS(theta) >= 2.0D0*bnd) THEN
  
  ind = 1
  CALL dcrec(REAL(t), AIMAG(t),u,v)
  c2 = -j * r * CMPLX(u,v)
  IF (is >= 0.0D0) THEN
    IF (is > 0.0D0) GO TO 10
    IF (REAL(s) <= 0.0D0) GO TO 10
  END IF
  c2 = -c2
END IF

!          SUMMATION OF THE SERIES S1 AND S2

10 sm = s * s
p = (a(2)*eta + a(3)) * s
p1 = ((a(4)*eta + a(5))*eta + a(6)) * sm
s1 = (1.0D0 + p) + p1
IF (ind /= 0) s2 = (1.0D0-p) + p1
sgn = 1.0D0
am = ar * ar
m = 4
l = 6

!          P = VALUE OF THE M-TH POLYNOMIAL

20 l = l + 1
alpha = a(l)
p = CMPLX(a(l),0.0D0)
DO  k = 2, m
  l = l + 1
  alpha = a(l) + aq * alpha
  p = a(l) + eta * p
END DO

!          ONLY THE S1 SUM IS FORMED WHEN IND = 0

sm = s * sm
p = p * sm
s1 = s1 + p
IF (ind /= 0) THEN
  sgn = -sgn
  s2 = s2 + sgn * p
END IF
am = ar * am
IF (1.0D0 + alpha*am /= 1.0D0) THEN
  m = m + 1
  IF (m <= 16) GO TO 20
END IF

!          FINAL ASSEMBLY

s1 = c * c1 * s1
IF (ind == 0) THEN
  w = s1
ELSE
  
  s2 = c * c2 * s2
  q = nu + z
  theta = ATAN2(AIMAG(q),REAL(q)) - phi
  IF (ABS(theta) <= bnd) THEN
    w = s1 + s2
  ELSE
    
    alpha = pi2
    IF (izn < 0.0D0) alpha = -alpha
    t = alpha * CMPLX(ABS(inu),-fn)
    alpha = EXP(REAL(t))
    u = AIMAG(t)
    r = CMPLX(COS(u),SIN(u))
    t = s1 - (alpha*r) * s1
    IF (x == 0.0D0 .AND. inu == 0.0D0) t = -t
    
    IF (y < 0.0D0) THEN
      IF (izn >= 0.0D0 .AND. theta <= SIGN(pi,theta)) s2 =  &
                                                      s2 * ( CONJG(r)/alpha)
      IF (x == 0.0D0) GO TO 40
      IF (izn >= 0.0D0) THEN
        IF (is < 0.0D0) GO TO 40
      END IF
    END IF
    
    w = s2 + t
    GO TO 50
    40 w = s2 - t
  END IF
END IF

50 IF (inu < 0.0D0) w = CONJG(w)
RETURN
END SUBROUTINE cbdb

!-----------------------------------------------------------------------------!


elemental SUBROUTINE dcrec(x, y, u, v)
!-----------------------------------------------------------------------
!             COMPLEX RECIPROCAL U + I*V = 1/(X + I*Y)
!-----------------------------------------------------------------------
REAL , INTENT(IN)   :: x, y
REAL , INTENT(OUT)  :: u, v

! Local variables
REAL   :: t, d

IF (ABS(x) <= ABS(y)) THEN
  t = x / y
  d = y + t * x
  u = t / d
  v = -1.0d0 / d
  RETURN
END IF
t = y / x
d = x + t * y
u = 1.0d0 / d
v = -t / d
RETURN
END SUBROUTINE dcrec

!-----------------------------------------------------------------------------!


end module Bessel_Functions
