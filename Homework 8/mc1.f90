  !  simple Monte Carlo integration of sin(x) 
  !
program mc1
  implicit none

  real(kind=8) :: x,iavg,i2avg,err,f
  integer,parameter :: n=1000000000
  real(kind=8),parameter :: exact=1.0d0-cos(1.0d0)
  integer :: i

  iavg=0.0d0
  i2avg=0.0d0

  do i=1,n
     call random_number(x)
     f=sin(x)

     iavg=iavg+f
     i2avg=i2avg+f*f

     if (mod(i,1000)==0) then
        print *,i,iavg/i,stderr(iavg,i2avg,i),abs(iavg/i-exact)
     endif
  enddo

contains

  ! calculate estimated standard error
  !
  ! xavg = sum of x_i
  ! x2avg = sum of x_i **2
  ! n = number of samples
  !
  function stderr(xavg,x2avg,n)
    real(kind=8) :: xavg,x2avg,stderr
    integer :: n
    real(kind=8) :: x,x2

    x=xavg/n
    x2=x2avg/n
    stderr=sqrt(1.0d0/(n-1)*(x2-x*x))
  end function stderr

end program mc1
