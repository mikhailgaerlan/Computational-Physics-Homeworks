! integrate.f90: Integrate exp(-x) using trapezoid, Simpson and Gauss rules 
!
! From: "A SURVEY OF COMPUTATIONAL PHYSICS" 
!	by RH Landau, MJ Paez, and CC BORDEIANU 
!	Copyright Princeton University Press, Princeton, 2008.
!	Electronic Materials copyright: R Landau, Oregon State Univ, 2008;
!	MJ Paez, Univ Antioquia, 2008; and CC BORDEIANU, Univ Bucharest, 2008.
!	Supported by the US National Science Foundation
!
! code cleaned up and rewritten in modern fortran style by RT Clay, 2013
!		
program integrate
  implicit none

  real(kind=4) :: r1, r2, r3
  real(kind=4) :: theo, vmin, vmax
  integer :: i

  !
  ! theoretical result, integration range
  !
  theo = 0.25*(2.0 + sin(2.0))  ! use built-in exp(x) function
  vmin=0.0
  vmax=1.0

  open(10, File='single/integsingle.dat', Status='Unknown')

  ! calculate integral using both methods for steps = 3..501
  do i= 3, 2001 , 2
     r1=trapez(i, vmin, vmax)
     r1=abs(r1-theo)
     r2=simpson(i,vmin, vmax)
     r2=abs(r2-theo)
     r3=quad(i,vmin, vmax)
     r3=abs(r3-theo)
     write(10,*) i, r1, r2, r3			
  end do
  close(10)

contains

!
! the function we want to integrate
!
  function f(x)
    real(kind=4) :: f, x

    f=(cos(x))**2

  end function f

!
! trapezoid rule
!
  function trapez(i, min, max)
    integer :: i, n				
    real(kind=4) :: interval, min, max, trapez, x
    trapez=0		 
    interval = ((max-min) / (i-1))

    ! sum the midpoints
    do	n=2, (i-1)						
       x = interval * (n-1)
       trapez = trapez + f(x)*interval
    end do

    ! add the endpoints	 
    trapez = trapez+0.5*(f(min)+f(max))*interval

  end function trapez

!
! Simpson's rule
!
  function simpson(i, min, max)
    integer :: i, n				
    real(kind=4) :: interval, min, max, simpson, x

    simpson=0.0		 
    interval = ((max-min) / (i-1))

    ! loop for odd points
    do	n=2, (i-1), 2						
       x = interval * (n-1)
       simpson = simpson + 4.0*f(x)
    end do
    ! loop for even points
    do	n=3, (i-1), 2						
       x = interval * (n-1)
       simpson = simpson + 2.0*f(x)
    end do
    ! add the endpoints	 
    simpson = simpson+f(min)+f(max)
    simpson=simpson*interval/3.0

  end function simpson

!
! Gauss' rule
!
  function quad(n, min, max)
    integer :: n
    real(kind=4),dimension(n) :: w,x  ! note use of automatic arrays
    real(kind=4) ::  min, max, quad

    integer :: i, job

    quad=0.0
    job=0.0
    call gauss(n, job, min, max, x, w)
    do	i=1, n
       quad=quad+f(x(i))*w(i)
    end do

  end function quad

!
!gauss.f90: Points and weights for Gaussian quadrature								 
!	rescale rescales the gauss-legendre grid points and weights
!
!	npts		 number of points
!	job =   0	 rescaling uniformly between (a,b)
!    	        1	 for integral (0,b) with 50% points inside (0, ab/(a+b))
!		2	 for integral (a,inf) with 50% inside (a,b+2a)
!		x, w		 output grid points and weights.
!
  subroutine gauss(npts,job,a,b,x,w) 
    integer,intent(in) ::npts,job
    real(kind=4),intent(in) :: a,b
    real(kind=4),intent(out) :: x(npts),w(npts)

    real(kind=4) :: xi,t,t1,pp,p1,p2,p3,aj
    real(kind=4),parameter :: pi=3.14159265358979323846264338328
    real(kind=4),parameter :: eps=3.0e-7
    real(kind=4),parameter :: zero=0.0, one=1.0, two=2.0
    real(kind=4),parameter :: half=0.5, quarter=0.25
    integer :: m,i,j 

    m=(npts+1)/2
    do	i=1,m
       t=cos(pi*(i-quarter)/(npts+half))
       do 
          p1=one
          p2=zero
          aj=zero
          do j=1,npts
             p3=p2
             p2=p1
             aj=aj+one
             p1=((two*aj-one)*t*p2-(aj-one)*p3)/aj
          end do
          pp=npts*(t*p1-p2)/(t*t-one)
          t1=t
          t=t1-p1/pp
          if (abs(t-t1)<eps) exit
       enddo
       x(i)=-t
       x(npts+1-i)=t
       w(i)=two/((one-t*t)*pp*pp)
       w(npts+1-i)=w(i)
    end do

    !
    ! rescale the grid points 
    !
    select case(job)		
    case (0)
       ! scale to (a,b) uniformly
       do i=1,npts
          x(i)=x(i)*(b-a)/two+(b+a)/two
          w(i)=w(i)*(b-a)/two
       end do
       
    case(1) 
       ! scale to (0,b) with 50% points inside (0,ab/(a+b))
       do i=1,npts
          xi=x(i)
          x(i)=a*b*(one+xi)/(b+a-(b-a)*xi)
          w(i)=w(i)*two*a*b*b/((b+a-(b-a)*xi)*(b+a-(b-a)*xi))
       end do
				 
    case(2) 
       ! scale to (a,inf) with 50% points inside (a,b+2a)
       do i=1,npts
          xi=x(i)
          x(i)=(b*xi+b+a+a)/(one-xi)
          w(i)=w(i)*two*(a+b)/((one-xi)*(one-xi))
       end do
    end select

  end subroutine gauss

 
end program integrate
