program mc1
  implicit none

  real(kind=8), dimension(10) :: x, y
  real(kind=8) :: iavg,isqavg,i2avg,i2sqavg,err,f,f2,rand,weight,int,int2
  integer,parameter :: n=20000000
  real(kind=8),parameter :: pi=4.0d0*atan(1.0d0), exact=(2.0d0*pi)**5
  integer :: i, j
  
  iavg=0.0d0
  isqavg=0.0d0
  i2sqavg=0.0d0
  i2sqavg=0.0d0
  
  do i=1,n
     do j=1,10
        call random_number(rand)
        x(j)=10.0d0*(rand-0.5d0)
        y(j)=-log(1.0d0-rand)
     enddo
     f=exp(-0.5d0*sum(x**2))
     weight=product(exp(-abs(y)))
     f2=exp(-0.5d0*sum(y**2))/weight
     
     iavg=iavg+f
     isqavg=isqavg+f*f
     i2avg=i2avg+f2
     i2sqavg=i2sqavg+f2*f2
     
     if (mod(i,1000)==0) then
        int = (1.0d10)*iavg/i
        int2 = ((2.0d0)**10)*i2avg/i
        print *,i,int,int2,stderr(iavg,isqavg,i),stderr(i2avg,i2sqavg,i),abs(int-exact)/exact,abs(int2-exact)/exact
     endif
  enddo

contains

  function stderr(xavg,x2avg,n)
    real(kind=8) :: xavg,x2avg,stderr
    integer :: n
    real(kind=8) :: x,x2

    x=xavg/n
    x2=x2avg/n
    stderr=sqrt(1.0d0/(n-1)*(x2-x*x))
  end function stderr

end program mc1
