program pendulum
  implicit none
  
  real(kind=8), parameter :: pi = 4.0_8*atan(1.0_8)
  
  integer :: n
  real(kind=8) :: t, theta  

  open(10, File='3/hw3-3.txt', Status='Unknown')
  
  theta = 0.0_8
  do while(theta < pi)
     t=0.0_8
     n=0
     call integ(theta, t, n)
     write(10,*) theta, nint(theta*180.0_8/pi), n, t, t/(2.0_8*pi)
     theta = theta + pi/45.0_8
  enddo
  
contains
  
  function f(x, theta)
    real(kind=8) :: f, x, theta
    
    f = 2.0_8/sqrt(1-(k(theta)**2)*(x**2))
  end function f

  function k(theta)
    real(kind=8) :: k, theta
    
    k = sin(theta/2.0_8)
  end function k

  subroutine integ(theta, t, n)
    integer, intent(out) :: n
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(out) :: t
    
    integer :: j
    real(kind=8) :: xj, wj, eps, oldt
    real(kind=8), parameter :: tol = 10.0_8**(-15)
    integer, parameter :: nmax = 1000
    
    do n=1,nmax
       oldt = t
       t=0.0_8
       do j=1,n
          xj=cos((pi*(real(j,8)-0.5_8))/real(n,8))
          wj=pi/real(n,8)
          t = t + wj*f(xj, theta)
       enddo
       if(abs(oldt-t)<tol)exit
    enddo
    if (n==nmax+1) then
       !print *, "Did not converge within ", nmax, " iterations."
       n = n - 1
    endif
    
  end subroutine integ
  
end program pendulum
