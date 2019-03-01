!
! solves ODE using 4th order Runge-Kutta
!
program rktest
  use deriv_mod
  use rk4_mod
  use rk2_mod
  implicit none
  
  real(kind=8), parameter :: pi = 4.0_8*atan(1.0_8)
  character(len=32) :: filename
  character(len=8) :: numstr
  real(kind=8), parameter :: maxt=5.0_8
  real(kind=8) :: dt
  integer, dimension(:), allocatable :: nstep
  
  real(kind=8), dimension(n) :: y
  real(kind=8) :: t
  integer :: i, j, s
  
  nstep = (/25, 30, 35, 100/)
  s = size(nstep)
  
  do j=1,s
     write(numstr, '(i0.2)') j
     filename = "data/data_rk4_"//trim(numstr)//".txt"
     open(10, file=filename, status="unknown")
     dt = maxt/nstep(j)
     write(10+j, *) maxt, nstep(j), dt, 0
     y(1)=1.0d0        ! initial conditions
     y(2)=0.0d0
     t=0.0d0           ! initialize independent variable
     
     do i=1,nstep(j)
        call rk4(t,dt,y)
        t=t+dt
        write(10,*) t, y(1), y(2), abs(y(1)-cos(2.0d0*pi*t))/abs(cos(2.0d0*pi*t)), &
             abs(y(2)+2.0d0*pi*sin(2.0d0*pi*t))/abs(2.0d0*pi*sin(2.0d0*pi*t))
     enddo
     
     close(10)
  end do
  
  do j=1,s
     write(numstr, '(i0.2)') j
     filename = "data/data_rk2_"//trim(numstr)//".txt"
     open(10, file=filename, status="unknown")
     dt = maxt/nstep(j)
     write(10+j, *) maxt, nstep(j), dt, 0
     y(1)=1.0d0        ! initial conditions
     y(2)=0.0d0
     t=0.0d0           ! initialize independent variable
     
     do i=1,nstep(j)
        call rk2(t,dt,y)
        t=t+dt
        write(10,*) t, y(1), y(2), abs(y(1)-cos(2.0d0*pi*t))/abs(cos(2.0d0*pi*t)), &
             abs(y(2)+2.0d0*pi*sin(2.0d0*pi*t))/abs(2.0d0*pi*sin(2.0d0*pi*t))
     enddo
     
     close(10)
  end do
  
end program rktest
