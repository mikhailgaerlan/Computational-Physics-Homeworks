program main
  implicit none
  
  real(kind=8), dimension(:,:), allocatable :: a
  real(kind=8), dimension(:), allocatable :: b
  real(kind=8), dimension(2,2) :: alpha, alu, ain
  real(kind=8), dimension(2) :: beta, x, ipiv, work
  integer :: i, info
  integer, parameter :: n = 20, m = 2
  real(kind=8), parameter :: sigma = 0.5d0
  
  allocate(a(n,m),b(n))
  
  open(11, file="data1.txt", status="old", action='read')
  
  do i = 1,n
     read(11,*), a(i,2), b(i)
     a(i,1) = 1.0d0/sigma
     a(i,2) = a(i,2)/sigma
     b(i) = b(i)/sigma
  enddo
  print *, " A                  b"
  print *, "------------------------"
  do i = 1, n
     print "(2f6.2,a,f6.2)", a(i,:), "       " ,b(i)
  enddo
  print *, ""
  
  call dgemm('T','N',m,m,n,1.0d0,a,n,a,n,1.0d0,alpha,m)
  print *, "alpha"
  print "(2f8.2)", alpha
  print *, ""
  
  call dgemv('T',n,m,1.0d0,a,n,b,1,1.0d0,beta,1)
  print *, "beta"
  print "(2f8.2)", beta
  print *, ""
  
  alu = alpha
  x = beta
  call dgesv(m,1,alu,m,ipiv,x,m,info)
  print *, "a"
  print "(2f8.4)", x
  print *, ""
  
  alu = alpha
  call dgetrf(m,m,alu,m,ipiv,info)
  ain = alu
  call dgetri(m,ain,m,ipiv,work,m,info)
  print *, "Uncertainties"
  print "(2f8.4)", sqrt(ain(1,1)), sqrt(ain(2,2))

end program main
