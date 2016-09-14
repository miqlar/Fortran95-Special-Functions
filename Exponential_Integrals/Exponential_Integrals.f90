!-----------------------------------------------------------------------!
!                          Exponential_Integrals                        !
!                                                                       !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Exponential_Integrals

use Incomplete_Gamma

implicit none

interface Ei
	module procedure Ei_R, Ei_C, En_IR, En_IC, En_RR, En_RC, En_CR, En_CC
end interface Ei

contains

!-----------------------------------------------------------------------!
!                           											!
!                                 Ei(z)                                 !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function Ei_R(x)
real, intent(in) :: x
real :: Ei_R, R, aux1, GA
integer :: i
real :: positive_infinity 

GA = 0.5772156649015328d0

positive_infinity=HUGE(positive_infinity)

if (x==0.0d0) then
	Ei_R=-positive_infinity
else if (x<40d0) then
	Ei_R=1.0d0
	aux1=1.0d0
	R=1.0d0
	do i=1,100
		R=R*i*x/((i+1.0d0)**2)
		Ei_R=Ei_R+R
	end do
	Ei_R=GA+LOG(x)+x*Ei_R
else
	Ei_R=1.0d0
	R=1.0d0
	do i=1,20
		R=R*i/x
		Ei_R=Ei_R+R
	end do
	Ei_R=EXP(x)/(x*Ei_R)
end if

end function Ei_R

!-----------------------------------------------------------------------!

elemental function Ei_C(z)
complex, intent(in) :: z
complex :: Ei_C

Ei_C = -Ei(1,-z)+0.5d0*LOG(z)-0.5d0*LOG(1.0d0/z)-LOG(-z)

end function Ei_C

!-----------------------------------------------------------------------!
!                           											!
!                                En(n,z)                                !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function En_IR(n,x)
real, intent(in) :: x
integer, intent(in) :: n
real :: En_IR, R, T0, GA, RP, PS, T, S, ENS, S0, positive_infinity
integer :: M, i, J, L, K, AllocateStatus
real, ALLOCATABLE:: aux_En_IR(:)

GA = 0.5772156649015328d0

positive_infinity=HUGE(positive_infinity)

if (n==1) then
	if (x==0.0d0) then
		En_IR=positive_infinity
	else if ( x<=1.0d0 ) then
		En_IR=1.0d0
		R=1.0d0
		do i=1,25
			R=-R*i*x/(i+1.0d0)**2
			En_IR=En_IR+R
		end do
		En_IR=-GA-log(x)+x*En_IR
	else
		M=20+INT(80.0/x)
		T0=0.0d0
		do i=M,1,-1
			T0=i/(1.0d0+i/(x+T0))
		end do
		T0=1.0d0/(x+T0)
		En_IR=EXP(-x)*T0
	end if
else 
    ALLOCATE(aux_En_IR(0:n), STAT = AllocateStatus)
   	if (x==0.0) then
		aux_En_IR(0)=(positive_infinity) 
		aux_En_IR(1)=(positive_infinity)
		do i=2,n
			aux_En_IR(i)=1.0d0/(i-1.0d0)
		end do
	else if (x<=1.0d0) then
		aux_En_IR(0)=(EXP(-x))/x
		do L=1,n
			RP=1.0d0
			do J=1,L-1
				RP=-RP*x/J
			end do
      PS =-0.5772156649015328d0
			do M=1,L-1
				PS=PS+1.0d0/M
			end do
			ENS=RP*(-LOG(x)+PS)
			S=0.0d0
			do M=0,20 
				if (M/=L-1) then
					R=1.0d0
					do J=1,M
						R=-R*x/J
					end do
					S=S+R/(M-L+1.0d0)
					S0=S
				end if
			end do
			aux_En_IR(L)=ENS-S
		end do	
	else
		aux_En_IR(0)=EXP(-x)/x
		M=15+INT(100.0/x)
		do L=1,n
			T0=0.0d0
			do K=M,1,-1
				T0=(L+K-1.0d0)/(1.0d0+K/(x+T0))
			end do
			T=1.0d0/(x+T0)
			aux_En_IR(L)=EXP(-x)*T
		end do
	end if
	En_IR=aux_En_IR(n)				
end if

end function En_IR

!-----------------------------------------------------------------------!

elemental function En_IC(n,z)
complex, intent(in) :: z
integer, intent(in) :: n
complex :: En_IC, CR, CT0, CT, sol
real    :: x, A0, EL, PI
integer :: K, positive_infinity, kerr

EL = 0.5772156649015328d0
PI = 3.141592653589793d0

positive_infinity=HUGE(positive_infinity)

x=real(z)
A0=abs(z)

if (n==1) then   ! E1(z)
	if (A0==0.0d0) then
		En_IC=CMPLX(positive_infinity,0.0d0) 
	else if (((A0<=10.0).OR.(x<0.0)).AND.A0<20.0) then
		En_IC=(1.0d0,0.0d0)
		CR=(1.0d0,0.0d0)
		do K=1,150
			CR=-CR*K*z/(K+1.0d0)**2
			En_IC=En_IC+CR
		end do
		En_IC=-EL-LOG(z)+z*En_IC
	else
		CT0=(0.0d0, 0.0d0)
		do K=120,1,-1
			CT0=K/(1.0d0+K/(Z+CT0))
		end do
		CT=1.0d0/(z+CT0)
		En_IC=EXP(-z)*CT
		if ((x<=0.0).AND.(AIMAG(z)==0.0d0)) then
			En_IC=En_IC-PI*(0.0d0,1.0d0)
		end if
	end if 
else	! EN(z)
  En_IC=(z**(n-1))*gamma(1-n, z)
end if

end function En_IC

!-----------------------------------------------------------------------!

elemental function En_RR(n,z)
real, intent(in) :: z
real, intent(in) :: n
complex :: c_En_RR, zz, nn
real::En_RR

zz=CMPLX(z,0)
nn=CMPLX(n,0)

c_En_RR=(zz**(nn-1))*gamma(1-nn, zz)

En_RR=real(c_En_RR)

end function En_RR

!-----------------------------------------------------------------------!

elemental function En_RC(n,z)
complex, intent(in) :: z
real, intent(in) :: n
complex :: En_RC, nn

nn=CMPLX(n,0)

En_RC=(z**(nn-1))*gamma(1-nn, z)

end function En_RC

!-----------------------------------------------------------------------!

elemental function En_CR(n,z)
real, intent(in) :: z
complex, intent(in) :: n
complex :: En_CR, zz

zz=CMPLX(z,0)

En_CR=(zz**(n-1))*gamma(1-n, zz)

end function En_CR

!-----------------------------------------------------------------------!

elemental function En_CC(n,z)
complex, intent(in) :: z
complex, intent(in) :: n
complex :: En_CC

En_CC=(z**(n-1))*gamma(1-n, z)

end function En_CC

!-----------------------------------------------------------------------!

end module Exponential_Integrals
