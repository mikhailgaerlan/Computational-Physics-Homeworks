! 
! Defines set of ODE's. In this example, N=2, simple harmonic motion
!
module deriv_mod
  implicit none

  private
  public:: deriv,n
  integer, parameter :: n=2

contains

  function deriv(x, y, i)
    real(kind=8), parameter :: pi = 4.0_8*atan(1.0_8)
    real(kind=8) :: deriv, m, k
    real(kind=8), intent(in) :: x     ! independent variable; not used here
    real(kind=8), intent(in) :: y(:)  ! dependent variable
    integer, intent(in) :: i          ! index of derivative needed
    
    m = 1.0_8
    k = 4.0_8*(pi**2)*m
    select case(i)
    case (1)
       deriv=y(2)
    case (2)
       deriv=-(k/m)*y(1)
    end select

  end function deriv

end module deriv_mod
