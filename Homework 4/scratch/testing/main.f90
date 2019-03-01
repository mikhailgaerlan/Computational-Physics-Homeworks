!
! solves ODE using 4th order Runge-Kutta
!
program rktest
  use deriv_mod
  use rk4_mod
  implicit none
  
  real(kind=8), parameter :: dt=0.01d0
  integer, parameter :: nstep=4000
  real(kind=8), dimension(3), parameter :: vi = (/0.1d0, 0.2d0, 0.25d0/)
  
  character(len=32) :: filename
  character(len=8) :: numstr, numstr2
  real(kind=8), dimension(n) :: y
  real(kind=8) :: t, b, y0, kinetic, potential, total
  integer :: i, m, j
  
  do j = 3, 3
     m = 1
     b = -1.0d0
     do while(b < 0)
        y0 = 5.0d0
        y(1) = b        ! initial conditions
        y(2) = -y0
        y(3) = 0.0d0
        y(4) = vi(j)
        t=0.0d0           ! initialize independent variable
        
        write(numstr, '(i0.3)') m
        write(numstr2, '(i0.1)') j
        filename = "data/data_"//trim(numstr2)//"_"//trim(numstr)//".txt"
        open(10, file = filename, status = 'unknown')
        
        do i=1,nstep
           call rk4(t,dt,y)
           t=t+dt
           kinetic = (1.0d0/2.0d0)*(0.5d0)*((y(3))**2)+(1.0d0/2.0d0)*(0.5d0)*((y(4))**2)
           potential = ((y(1))**2)*((y(2))**2)*exp(-((y(1))**2+(y(2))**2))
           total = kinetic + potential
           write(10,*) t, y(1), y(2), y(3), y(4), kinetic, potential, total, atan2(y(4),y(3)), potential/kinetic
        enddo
        
        b = b + 0.01d0
        m = m + 1
     enddo
  enddo
end program rktest
