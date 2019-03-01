!
! solves ODE using 4th order Runge-Kutta
!
program rktest
  use deriv_mod
  use rk4_mod
  implicit none

  real(kind=8),parameter :: dt=0.01d0
  integer,parameter :: nstep=1000
  
  real(kind=8),dimension(n) :: y
  real(kind=8) :: t
  integer :: i

  y(1)=1.0d0        ! initial conditions
  y(2)=0.0d0
  t=0.0d0           ! initialize independent variable
  
  do i=1,nstep
     call rk4(t,dt,y)
     t=t+dt
     print *,t,y(1),y(2)
  enddo

end program rktest
