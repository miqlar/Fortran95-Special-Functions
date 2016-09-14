!-----------------------------------------------------------------------!
!                                                                       !
!                            Error_Function                             !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Error_Function

implicit none

interface erf
	module procedure Error_real, Error_complex
end interface erf

interface erfc
	module procedure c_Error_real, c_Error_complex
end interface erfc

contains

!-----------------------------------------------------------------------!

elemental function Error_real(x)
real, intent(in) :: x
real :: EPS, pi, ER, x2, R, c0, ERR
real :: Error_real
integer :: k

pi = 4d0*atan(1d0)
x2=x*x
ER=1.0d0
R=1.0d0

if (ABS(x)<3.5d0) then
	do k=1,50 
		R=R*x2/(K+0.5d0)
		ER=ER+R
	end do
	C0=2.0d0/SQRT(pi)*x*EXP(-x2)
	Error_real=C0*ER
else
	do k=1,12
		R=-R*(K-0.5d0)/x2
		ER=ER+R
	end do
	C0=EXP(-x2)/(ABS(x)*SQRT(pi))
	Error_real=1.0d0-C0*ER
	if (x<0.0) then
		Error_real=-Error_real
	end if
end if 

end function Error_real
!-----------------------------------------------------------------------!

elemental function Error_complex(z)
complex, intent(in) :: z
integer :: k
real :: A0, pi
complex :: z1, CS, CR, CL, C0
complex :: Error_complex

A0=abs(z)
C0=EXP(-z*z)
pi = 4d0*atan(1d0)
z1=z

if (real(z)<0.0d0) then
	z1=-z
end if 
if (A0<=5.8d0) then
	CS=z1
	CR=z1
	do k=1,120
		CR=CR*z1*z1/(k+0.5d0)
		CS=CS+CR
	end do
	Error_complex=2.0d0*C0*CS/SQRT(pi)
else
	CL=1.0d0/z1
	CR=CL
	do k=1,13
		CR=-CR*(k-0.5d0)/(z1*z1)
		CL=CL+CR
	end do
	Error_complex=1.0d0-C0*CL/SQRT(pi)
end if
if (real(z)<0.0d0) then
	Error_complex=-Error_complex
end if 

end function Error_complex

!-----------------------------------------------------------------------!

elemental function c_Error_real(x)
real, intent(in) :: x
real :: c_Error_real

c_Error_real = 1-erf(x)

end function c_Error_real

!-----------------------------------------------------------------------!

elemental function c_Error_complex(z)
complex, intent(in) :: z
complex :: c_Error_complex

c_Error_complex = 1-erf(z)

end function c_Error_complex

!-----------------------------------------------------------------------!

end module Error_Function