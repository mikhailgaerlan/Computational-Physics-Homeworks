!
! solves ODE using 4th order Runge-Kutta
!
program rktest
  use deriv_mod
  use rk4_mod
  use rk2_mod
  implicit none
  
  real(kind=8), parameter :: dt = pi/100.0d0
  integer, parameter :: nstep = 30025
  
  real(kind=8), dimension(n) :: y
  real(kind=8) :: t
  integer :: i
  
  open(10, file='data/data.txt',status='unknown')
  f = 0.0d0
  do while(f < 2.25)
     f = f + 0.0001
     y(1) = 0.0d0
     y(2) = 1.0d0
     t = 0.0d0
     do i=1,nstep
        call rk4(t,dt,y)
        t=i*dt
        if (y(1)<-pi) y(1)=y(1)+2*pi
        if (y(1)>pi) y(1)=y(1)-2*pi
     enddo
     write(10,*) f, abs(y(2))
  enddo
  close(10)
  
end program rktest
