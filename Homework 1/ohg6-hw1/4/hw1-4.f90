program homework
  implicit none

  integer :: i, n

  real(kind=8) :: x1, term1, sum1, x2, term2, sum2
  real(kind=8), parameter :: eps=1.0e-10_8

  print *, " x    n   (a)                   (b)"
  print *, "------------------------------------------------"
  do i=1,10
     x1=real(i,8)*10.0_8
     x2=-x1
     n=0
     term1=1.0_8
     term2=term1
     sum1=term1
     sum2=term2
     do
        n=n+1
        term1=-term1*x1/n
        term2=-term2*x2/n
        if (abs(term1)<eps .or. abs(term2)<eps) exit
        sum1=sum1+term1
        sum2=sum2+term2
     enddo
     print "(i3,a,i3,es20.10,a,es20.10)", nint(x1),"  ", n, sum1," ", 1.0_8/sum2
  enddo
  
end program homework
