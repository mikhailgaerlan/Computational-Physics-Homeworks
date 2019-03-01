program main
  implicit none
  
  real(kind=8), dimension(:,:), allocatable :: a
  real(kind=8), dimension(:), allocatable :: b
  real(kind=8), dimension(3,3) :: alpha, alu, ain
  real(kind=8), dimension(3) :: beta, x, ipiv, work
  real(kind=8) :: theta, rho
  integer :: i, info
  real(kind=8), parameter :: pi = 4.0d0*atan(1.0d0)
  integer, parameter :: n = 10, m = 3
  
  allocate(a(n,m),b(n))
  
  open(11, file="data2.txt", status="old", action='read')
  open(12, file="dataout.txt", status="unknown", action='write')
  
  do i = 1,n
     b(i) = 1.0d0
     read(11,*), a(i,1), a(i,2)
     a(i,3) = a(i,1)*a(i,2)
     a(i,1) = (a(i,1))**2
     a(i,2) = (a(i,2))**2
  enddo
  print *, " A                            b"
  print *, "----------------------------------"
  do i = 1,n
     print "(3f8.4,a,f4.2)", a(i,1:3), "       " ,b(i)
  enddo
  print *, ""
  
  call dgemm('T','N',m,m,n,1.0d0,a,n,a,n,1.0d0,alpha,m)
  print *, "alpha"
  print "(3f9.4)", alpha
  print *, ""
  
  call dgemv('T',n,m,1.0d0,a,n,b,1,1.0d0,beta,1)
  print *, "beta"
  print "(3f9.4)", beta
  print *, ""
  
  alu = alpha
  x = beta
  call dgesv(m,1,alu,m,ipiv,x,m,info)
  print *, "a"
  print "(3f9.4)", x
  print *, ""

  do i=1,1000
     theta = 2*pi*i/1000.0d0
     rho = 1.0d0/sqrt(x(1)*(cos(theta))**2+x(2)*(sin(theta))**2+x(3)*(sin(theta)*cos(theta)))
     write(12,*) rho*cos(theta), rho*sin(theta)
  enddo
  
end program main
