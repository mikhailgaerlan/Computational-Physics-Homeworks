!  From: "A SURVEY OF COMPUTATIONAL PHYSICS" 
!  by RH Landau, MJ Paez, and CC BORDEIANU 
!  Copyright Princeton University Press, Princeton, 2008.
!  Electronic Materials copyright: R Landau, Oregon State Univ, 2008;
!  MJ Paez, Univ Antioquia, 2008; and CC BORDEIANU, Univ Bucharest, 2008.
!  Support by National Science Foundation  		 
! 		
!   rk45.f90: ODE solver via variable step size rk, Tol = error
!
! cleaned up and rewritten in modern F90 style  by RTC 09/2013
program rk45 
  implicit none

  real(kind=8) :: h, t, s, s1, hmin, hmax
  real(kind=8),parameter :: Tol = 2.0d-7
  real(kind=8),parameter :: Tmin = 0.0d0
  real(kind=8),parameter :: Tmax = 10.0d0
  real(kind=8),dimension(2) :: w, y, FReturn, ydumb, k1, k2, k3, k4, k5, k6, err
  integer :: i
  integer,parameter :: Ntimes = 10  

  ! output to file	
  Open(10, FILE = 'rk45.dat', Status = 'Unknown')
  
  ! initialize
  y(1) = 3.0 ; y(2) = - 5.0	

  ! tentative number of steps			                              
  h = (Tmax - Tmin) / Ntimes 
    
  ! minimum and maximum step size
  hmin = h/64    
  hmax = h*64										      
  t = Tmin

  do while (t < Tmax)
     write(10, *) t, y(1)	                                    
     If ( (t + h)  >  Tmax ) then 
        h = Tmax - t	! the last step
     endif
     ! evaluate both RHSs and Return in F
     call f(t, y, FReturn)		            
     do i = 1, 2 
        k1(i) = h*FReturn(i)		
        ydumb(i) = y(i) + k1(i)/4
     end do
     call f(t + h/4, ydumb, FReturn)
     do i = 1, 2 
        k2(i) = h*FReturn(i)
        ydumb(i) = y(i) + 3.0d0*k1(i)/32 + 9.0d0*k2(i)/32
     end do
     call f(t + 3.0d0*h/8, ydumb, FReturn)
     do i = 1, 2 
        k3(i) = h*FReturn(i)	
        ydumb(i) = y(i) + 1932.0d0*k1(i)/2197.0d0 - 7200.0d0*k2(i)/2197.0d0 &
             + 7296.0d0*k3(i)/2197.0d0
     end do
     call f(t + 12*h/13.0d0, ydumb, FReturn)
     do i = 1, 2
        k4(i) = h*FReturn(i) 
        ydumb(i) = y(i) + 439.0d0*k1(i)/216.0d0-8.0d0*k2(i) &
             + 3680.0d0*k3(i)/513.0d0-845.0d0*k4(i)/4104.0d0
     end do
     call f(t + h, ydumb, FReturn)
     do i = 1, 2 
        k5(i) = h*FReturn(i)
        ydumb(i) = y(i) - 8*k1(i)/27.0d0 + 2*k2(i) - 3544.0d0*k3(i)/2565.0d0 &
             + 1859.0d0*k4(i)/4104.0d0 - 11.0d0*k5(i)/40.0d0
     end do
     call f(t + h/2, ydumb, FReturn)
     do i = 1, 2 
        k6(i) = h*FReturn(i) 
        err(i) = abs( k1(i)/360.0d0 - 128*k3(i)/4275.0d0 - 2197.0d0*k4(i)/75240.0d0 &
             + k5(i)/50.0d0 + 2*k6(i)/55.0d0 )
     end do
     If ((err(1) < Tol).or.(err(2) < Tol).or.(h <= 2*hmin))	then 
        ! accept approximation
        do i = 1, 2                                    
           y(i) = y(i) + 25.0d0*k1(i)/216.0d0 + 1408.0d0*k3(i)/2565.0d0 &
                + 2197.0d0*k4(i)/4104.0d0 - k5(i)/5.0d0
        end do
        t = t + h
     endif
     If (( err(1) == 0).or.(err(2) == 0)) then
        s = 0.0d0! trap division by 0
     else 
        s = 0.84d0*Tol*h/err(1)**0.25d0                  ! step size scalar
     endif
     If ( (s  <  0.75d0).and. (h  >  2*hmin) )then
        h = 	h/2.0d0 ! reduce step
     else If ( (s  >  1.5d0).and.(2* h  <  hmax) )then
        h = h*2.0d0	! increase step
     endif
  end do
  close(10)
  stop'Data stored in rk45.dat'

contains
        
  ! PLACE YOUR FUNCTION HERE
  subroutine f(t,y,FReturn) 
    real(kind=8) :: t, y(2), FReturn(2)

    FReturn(1) = y(2)! RHS of first equation
    FReturn(2) = - 100.0d0*y(1) - 2.0d0*y(2) + 10.0d0*sin(3.0d0*t)! RHS of 2nd equation

  end subroutine f

end program rk45

