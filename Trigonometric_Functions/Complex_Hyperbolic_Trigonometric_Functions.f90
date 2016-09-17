!-----------------------------------------------------------------------!
!              Complex_Hyperbolic_Trigonometric_Functions               !
!                                                                       !
!                                                                       !
!--------------------------------------------Miquel Larsson---2015------!

module Complex_Hyperbolic_Trigonometric_Functions
implicit none

intrinsic sinh
intrinsic cosh
intrinsic tanh

interface sinh
  module procedure csinhc, csinhr
end interface sinh

interface cosh
  module procedure ccoshc, ccoshr
end interface cosh

interface tanh
  module procedure ctanhc, ctanhr
end interface tanh


contains

!-----------------------------------------------------------------------!
!                           											!
!                              sinh(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function csinhc(z)
complex, intent(in) :: z
complex             :: csinhc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

csinhc = sinh(x)*cos(y)+i*cosh(x)*sin(y)

end function csinhc

!-----------------------------------------------------------------------!

elemental function csinhr(z)
real, intent(in) :: z
complex          :: csinhr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

!csinhr = (sinh(x)*cos(y))+(i*cosh(x)*sin(y)) RECURSION ISSUE
csinhr=-i*(sin(i*x))

end function csinhr

!-----------------------------------------------------------------------!
!                           											!
!                              cosh(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function ccoshc(z)
complex, intent(in) :: z
complex             :: ccoshc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

ccoshc=cosh(x)*cos(y)+i*sinh(x)*sin(y)

end function ccoshc

!-----------------------------------------------------------------------!

elemental function ccoshr(z)
real, intent(in) :: z
complex          :: ccoshr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

ccoshr=cosh(x)*cos(y)+i*sinh(x)*sin(y)

end function ccoshr

!-----------------------------------------------------------------------!
!                           											!
!                              tanh(z)                                  !
!                                                                       !
!-----------------------------------------------------------------------!

elemental function ctanhc(z)
complex, intent(in) :: z
complex             :: ctanhc

complex :: i
real    :: x, y

i = (0d0,1d0)
x = real(z)
y = aimag(z)

ctanhc=(tanh(2d0*x)/(1d0+cos(2d0*y)/cosh(2d0*x)))+(i*sin(2d0*y)/(cosh(2d0*x)+cos(2d0*y)))

end function ctanhc

!-----------------------------------------------------------------------!

elemental function ctanhr(z)
real, intent(in) :: z
complex          :: ctanhr

complex :: i
real    :: x, y

i = (0d0,1d0)
x = z
y = 0d0

ctanhr=(tanh(2d0*x)/(1d0+cos(2d0*y)/cosh(2d0*x)))+(i*sin(2d0*y)/(cosh(2d0*x)+cos(2d0*y)))

end function ctanhr

!-----------------------------------------------------------------------!

end module Complex_Hyperbolic_Trigonometric_Functions
