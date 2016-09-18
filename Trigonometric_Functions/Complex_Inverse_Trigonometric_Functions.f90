!-----------------------------------------------------------------------!
!                 Complex_Inverse_Trigonometric_Functions               !
!                                                                       !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Complex_Inverse_Trigonometric_Functions
implicit none

real :: pi = 3.1415926535897932

intrinsic asin
intrinsic acos
intrinsic atan

interface asin
  module procedure casinc, casinr
end interface asin

interface acos
  module procedure cacosc, cacosr
end interface acos

interface atan
  module procedure catanc, catanr
end interface atan


contains

!-----------------------------------------------------------------------!
!                           											!
!                              asin(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function casinc(z)
complex, intent(in) :: z
complex             :: casinc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

casinc =     0.5d0 * sign(1d0,x) * acos( sqrt((abs(z)**2 - 1d0)**2 + 4d0 * y**2) - abs(z)**2) + &
         i * 0.5d0 * sign(1d0,y) * acosh(sqrt((abs(z)**2 - 1d0)**2 + 4d0 * y**2) + abs(z)**2)

end function casinc

!-----------------------------------------------------------------------!

elemental function casinr(z)
real, intent(in) :: z
complex          :: casinr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

casinr =     0.5d0 * sign(1d0,x) * acos( sqrt((abs(z)**2 - 1d0)**2 + 4d0 * y**2) - abs(z)**2) + &
         i * 0.5d0 * sign(1d0,y) * acosh(sqrt((abs(z)**2 - 1d0)**2 + 4d0 * y**2) + abs(z)**2)

end function casinr

!-----------------------------------------------------------------------!
!                           											!
!                              acos(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function cacosc(z)
complex, intent(in) :: z
complex             :: cacosc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

cacosc = (pi/2d0) - asin(z)

end function cacosc

!-----------------------------------------------------------------------!

elemental function cacosr(z)
real, intent(in) :: z
complex             :: cacosr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

cacosr = (pi/2d0) - asin(z)

end function cacosr

!-----------------------------------------------------------------------!
!                           											!
!                              atan(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function catanc(z)
complex, intent(in) :: z
complex             :: catanc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

if (x/=0) then
	catanc = 0.5d0*atan2(2d0*x,1-abs(z)**2)+i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
else if (abs(y)>1) then
	catanc = 0.5d0*pi*sign(1d0,y)+i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
else 
	catanc = i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
end if

end function catanc

!-----------------------------------------------------------------------!

elemental function catanr(z)
real, intent(in) :: z
complex          :: catanr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

if (x/=0) then
	catanr = 0.5d0*atan2(2d0*x,1d0-abs(z)**2)+i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
else if (abs(y)>1) then
	catanr = 0.5d0*pi*sign(1d0,y)+i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
else 
	catanr = i*0.5d0*atanh(2d0*y/(1d0+abs(z)**2))
end if

end function catanr

!-----------------------------------------------------------------------!

end module Complex_Inverse_Trigonometric_Functions
