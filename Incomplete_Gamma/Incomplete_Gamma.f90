!-----------------------------------------------------------------------!
!                           Incomplete_Gamma                            !
!                                                                       !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Incomplete_Gamma
implicit none

integer, parameter  :: dp = SELECTED_REAL_KIND(12, 60)

interface Gamma
  module procedure Inc_Gamma_IC, Inc_Gamma_CC, Inc_gamma_RC
end interface Gamma

contains

!----------------------------------------------------------------------!
!                               Inc_Gamma_IC                           !
!                                                                      !
!    Incomplete Gamma function for an integer and a complex argument   !
!                                                                      ! 
!                                                                      ! 
!--------------------------------------------Miquel Larsson-2015.03.18-!

elemental function Inc_Gamma_IC(n,z)
complex, intent(in) :: z
integer, intent(in) :: n
complex :: Inc_Gamma_IC, aux

aux=CMPLX(n,0)
Inc_Gamma_IC=cdig(aux,z)

end function Inc_Gamma_IC

!----------------------------------------------------------------------!
!                               Inc_Gamma_CC                           !
!                                                                      !
!          Incomplete Gamma function for two complex arguments         !
!                                                                      ! 
!                                                                      ! 
!--------------------------------------------Miquel Larsson-2015.03.18-!

elemental function Inc_Gamma_CC(n,z)
complex, intent(in) :: z
complex, intent(in) :: n
complex :: Inc_Gamma_CC

Inc_Gamma_CC=cdig(n,z)

end function Inc_Gamma_CC

!----------------------------------------------------------------------!
!                               Inc_Gamma_RC                           !
!                                                                      !
!     Incomplete Gamma function for a real and a complex argument      !
!                                                                      ! 
!                                                                      ! 
!--------------------------------------------Miquel Larsson-2015.03.18-!

elemental function Inc_Gamma_RC(n,z)
complex, intent(in) :: z
real, intent(in) :: n
complex :: Inc_Gamma_RC, aux

aux=CMPLX(n,0)
Inc_Gamma_RC=cdig(aux,z)

end function Inc_Gamma_RC

!-----------------------------------------------------------------------!
!                           											!
!                Functions that compute the incomplete                  !
!           gamma function for complex argument and order.              !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function cdig(alpha, x) RESULT(fn_val)
 
complex , intent(IN) :: alpha
complex , intent(IN) :: x
complex              :: fn_val

complex              :: p, q
integer              :: i, ilim
real , parameter     :: zero = 0.0_dp, xlim = 1.0_dp
complex , parameter  :: cone = (1.0_dp, 0.0_dp)
real , parameter     :: re = (0.36787944117144232_dp, 0.0_dp)
integer, parameter   :: ibuf = 34

! --- If x is near the negative real axis, then shift to x=1.
IF (dnrm(x) < xlim .OR. real(x, KIND=dp) < zero .AND. ABS(AIMAG(x)) < xlim) THEN
  fn_val = re / cdh(alpha, cone)
  ilim = real(x/re, KIND=dp)
  DO  i = 0, ibuf - ilim
    CALL term(alpha, x, i, p, q)
    fn_val = fn_val + p * q
  END DO
ELSE
  fn_val = EXP(-x + alpha*LOG(x)) / cdh(alpha, x)
END IF
RETURN
END function cdig

!-----------------------------------------------------------------------!

elemental function cdh(alpha, x) RESULT(fn_val)
! --- Written By Eric Kostlan & Dmitry Gokhman
! --- March  1986

complex , intent(IN)  :: alpha
complex , intent(IN)  :: x
complex               :: fn_val

complex   :: term, sum, cn, alpha1
integer       :: i, n
real , parameter  :: one = 1.0_dp

! --- If Re(alpha-x) is too big, shift alpha.
n = real(alpha-x, KIND=dp)
IF (n > 0) THEN
  cn = n
  alpha1 = alpha - cn
  term = one / x
  sum = term
  DO  i = 1, n - 1
    cn = n - i
    term = term * (alpha1 + cn) / x
    sum = term + sum
  END DO
  sum = sum + term * alpha1 / cdhs(alpha1, x)
  fn_val = one / sum
ELSE
  fn_val = cdhs(alpha, x)
END IF
RETURN
END function cdh

!-----------------------------------------------------------------------!


elemental function cdhs(alpha, x) RESULT(fn_val)
! --- Written By Eric Kostlan & Dmitry Gokhman
! --- March  1986

complex , intent(IN)  :: alpha
complex , intent(IN)  :: x
complex               :: fn_val

complex   :: p0, q0, p1, q1, r0, r1, ci, factor
integer       :: i
real , parameter  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp
real , parameter  :: tol1 = 1.0D+10, tol2 = 1.0D-10, error = 5.D-18
integer, parameter    :: ilim = 100000

q0 = one
q1 = one
p0 = x
p1 = x + one - alpha
DO  i = 1, ilim
  ci = i
  IF (p0 /= zero .AND. q0 /= zero .AND. q1 /= zero) THEN
    r0 = p0 / q0
    r1 = p1 / q1
    IF (dnrm(r0-r1) <= dnrm(r1)*error) THEN
      fn_val = r1
      RETURN
    END IF
! --------- Occasionally renormalize the sequences to avoid over(under)flow.
    IF (dnrm(p0) > tol1 .OR. dnrm(p0) < tol2 .OR. dnrm(q0) > tol1  &
          .OR. dnrm(q0) < tol2) THEN
      factor = p0 * q0
      p0 = p0 / factor
      q0 = q0 / factor
      p1 = p1 / factor
      q1 = q1 / factor
    END IF
  END IF
  p0 = x * p1 + ci * p0
  q0 = x * q1 + ci * q0
  p1 = p0 + (ci+one-alpha) * p1
  q1 = q0 + (ci+one-alpha) * q1
END DO
! --- If the peripheral routines are written correctly,
! --- the following four statements should never be executed.
!WRITE(*, *) 'cdhs:  *** Warning: i >', ilim
!WRITE(*, *) 'cdhs:  *** r0,r1= ', r0, r1
fn_val = half * (r0+r1)
RETURN
END function cdhs

!-----------------------------------------------------------------------!

elemental subroutine term(alpha, x, i, p, q)
! --- Calculate p*q = -1**i(1 - x**(alpha+i))/(alpha+i)i ! carefully.

complex , intent(IN)   :: alpha
complex , intent(IN)   :: x
integer, intent(IN)    :: i
complex , intent(OUT)  :: p
complex , intent(OUT)  :: q

complex   :: ci, cdlx, alphai
real , parameter  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
real , parameter  :: tol = 3.0D-7, xlim = 39.0_dp

IF (i == 0) q = one
ci = i
alphai = alpha + ci
IF (x == zero) THEN
  p = one / alphai
  IF (i /= 0) q = -q / ci
  RETURN
END IF
cdlx = LOG(x)

! --- If (1 - x**alphai) = -x**alphai on the computer,
! --- then change the inductive scheme to avoid overflow.
IF (real(alphai*cdlx, KIND=dp) > xlim .AND. i /= 0) THEN
  p = p * (alphai - one) / alphai
  q = -q * x / ci
  RETURN
END IF
IF (dnrm(alphai) > tol) THEN
  p = (one - x**alphai) / alphai
ELSE
  p = -cdlx * (one + cdlx*alphai/two)
END IF
IF (i /= 0) q = -q / ci
RETURN
END subroutine term

!-----------------------------------------------------------------------!

elemental function dnrm(z) RESULT(fn_val)

complex , intent(IN)  :: z
real                  :: fn_val

fn_val = ABS(real(z, KIND=dp)) + ABS(AIMAG(z))
RETURN
END function dnrm

end module Incomplete_Gamma
