! 
! Defines set of ODE's. In this example, N=2, simple harmonic motion
!
module deriv_mod
  implicit none

  private
  public:: deriv,n
  integer,parameter :: n=4

contains

  function deriv(x, y, i)
    real(kind=8) :: deriv, m
    real(kind=8),intent(in) :: x     ! independent variable; not used here
    real(kind=8),intent(in) :: y(:)  ! dependent variable
    integer,intent(in) :: i          ! index of derivative needed
    
    m = 0.5d0
    select case(i)
    case (1)
       deriv = y(3)
    case (2)
       deriv = y(4)
    case(3)
       deriv = (-1.0d0/m)*2.0d0*(y(2)**2)*y(1)*(1.0d0-y(1)**2)*exp(-(y(1)**2+y(2)**2))
    case(4)
       deriv = (-1.0d0/m)*2.0d0*(y(1)**2)*y(2)*(1.0d0-y(2)**2)*exp(-(y(1)**2+y(2)**2))
    end select

  end function deriv

end module deriv_mod
