! 
! Defines set of ODE's. In this example, N=2, simple harmonic motion
!
module deriv_mod
  implicit none

  private
  public:: deriv,n,f,alpha,pi,omega0,omega
  real(kind=8) :: f
  integer, parameter :: n=2
  real(kind=8), parameter :: pi = 4.0_8*atan(1.0_8), alpha = 0.1d0, omega0 = 1.0d0, omega = 2.0d0

contains

  function deriv(x, y, i)
    real(kind=8) :: deriv
    real(kind=8), intent(in) :: x     ! independent variable; not used here
    real(kind=8), intent(in) :: y(:)  ! dependent variable
    integer, intent(in) :: i          ! index of derivative needed
    
    select case(i)
    case (1)
       deriv=y(2)
    case (2)
       deriv=-alpha*y(2)-(omega0**2+f*cos(omega*x))*sin(y(1))
    end select

  end function deriv

end module deriv_mod
