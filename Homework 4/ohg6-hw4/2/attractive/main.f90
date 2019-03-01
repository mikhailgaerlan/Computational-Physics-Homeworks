!
! solves ODE using 4th order Runge-Kutta
!
program rktest
  use deriv_mod
  use rk4_mod
  implicit none
  
  real(kind=8), parameter :: dt=0.01d0
  integer :: nstep
  real(kind=8), dimension(3), parameter :: vi = (/0.5d0,1.0d0,1.5d0/)
  
  character(len=32) :: filename
  character(len=8) :: numstr, numstr2
  real(kind=8), dimension(n) :: y
  real(kind=8) :: t, b, y0, kinold, potold, kinetic, potential, total, angle, db
  integer :: i, m, j
  
  do j = 1, 3
     nstep = 3000
     write(numstr2, '(i0.1)') j
     m = 1
     b = -1.0d0
     do while(b < 0.0d0)
        y0 = 5.0d0
        y(1) = b        ! initial conditions
        y(2) = -y0
        y(3) = 0.0d0
        y(4) = vi(j)
        t=0.0d0           ! initialize independent variable
        
        write(numstr, '(i0.3)') m
        filename = "attractive_data/data_"//trim(numstr2)//"_"//trim(numstr)//".txt"
        
        open(10, file = filename, status = 'unknown')
        
        do i=1,nstep
           call rk4(t,dt,y)
           t=t+dt
           write(10,*) t, y(1), y(2), y(3), y(4)
        enddo
        
        b = b + 0.1d0
        m = m + 1
        close(10)
     enddo
     
     nstep = 10000
     b = -1.0d0
     db = 1.0d0/real(300,8)
     
     filename = "attractive_data/data_scatt_"//trim(numstr2)//".txt"
     open(10, file = filename, status = 'unknown')
     do while(b < 0.0d0)
        y0 = 5.0d0
        y(1) = b        ! initial conditions
        y(2) = -y0
        y(3) = 0.0d0
        y(4) = vi(j)
        t=0.0d0           ! initialize independent variable
        
        open(10, file = filename, status = 'unknown')
        
        kinold = (1.0d0/2.0d0)*(0.5d0)*((y(3))**2)+(1.0d0/2.0d0)*(0.5d0)*((y(4))**2)
        potold = ((y(1))**2)*((y(2))**2)*exp(-((y(1))**2+(y(2))**2))
        do i=1,nstep
           call rk4(t,dt,y)
           t=t+dt
        enddo
        kinetic = (1.0d0/2.0d0)*(0.5d0)*((y(3))**2)+(1.0d0/2.0d0)*(0.5d0)*((y(4))**2)
        potential = ((y(1))**2)*((y(2))**2)*exp(-((y(1))**2+(y(2))**2))
        angle = atan2(y(4),y(3))
        if (angle < 0) angle = angle + 8.0d0*atan(1.0d0)
        write(10, *) b, angle, kinold, potold, potold/kinold, kinetic, potential, potential/kinetic
        
        b = b + db
     enddo
     close(10)
  enddo
end program rktest
