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
  
  integer :: i
  
  real(kind=8), parameter :: pi=4.0_8*atan(1.0_8)
  real(kind=8) :: x
  !real(kind=8), dimension(4) :: z
  real(kind=8) :: z
  real(kind=8) :: cu, su
  
  open(10, File='1/hw3-1.txt', Status='Unknown')
  
  x = -5.0_8
  !z = (/ 0.5_8, 1.0_8, 2.0_8, 3.0_8 /)
  z = 0.5_8
  
  do while(z .le. 3)
     x=-5.0_8
     do while(x .le. 5)
        cu=quadc(100,0.0_8,u(x,z))
        su=quads(100,0.0_8,u(x,z))
        
        write(10, *) x, z, ((2*cu+1)**2+(2*su+1)**2)/8.0_8
        x=x+0.015625
     enddo
     z=z+0.015625
  enddo
  
contains
  
  function u(x,z)
    real(kind=8) :: u, x, z
    
    u = x*sqrt(2.0_8/z)
    
  end function u
  
!===================================
! the function we want to integrate
!===================================
  function c(x)
    real(kind=8) :: c, x

    c=cos(pi*(x**2)/2.0_8)

  end function c

  function s(x)
    real(kind=8) :: s, x

    s=sin(pi*(x**2)/2.0_8)

  end function s

!=============
! Gauss' rule
!=============
  function quadc(n, min, max)
    integer :: n
    real(kind=8),dimension(n) :: w,x  ! note use of automatic arrays
    real(kind=8) ::  min, max, quadc

    integer :: i, job

    quadc=0.0_8
    job=0.0_8
    call gauss(n, job, min, max, x, w)
    do	i=1, n
       quadc=quadc+c(x(i))*w(i)
    end do

  end function quadc

  function quads(n, min, max)
    integer :: n
    real(kind=8), dimension(n) :: w,x  ! note use of automatic arrays
    real(kind=8) ::  min, max, quads

    integer :: i, job

    quads=0.0_8
    job=0.0_8
    call gauss(n, job, min, max, x, w)
    do	i=1, n
       quads=quads+s(x(i))*w(i)
    end do

  end function quads

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
    real(kind=8),intent(in) :: a,b
    real(kind=8),intent(out) :: x(npts),w(npts)

    real(kind=8) :: xi,t,t1,pp,p1,p2,p3,aj
    real(kind=8),parameter :: pi=3.14159265358979323846264338328_8
    real(kind=8),parameter :: eps=3.0e-16_8
    real(kind=8),parameter :: zero=0.0_8, one=1.0_8, two=2.0_8
    real(kind=8),parameter :: half=0.5_8, quarter=0.25_8
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
