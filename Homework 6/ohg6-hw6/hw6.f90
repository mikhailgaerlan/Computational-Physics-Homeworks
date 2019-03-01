program main
  implicit none

  integer :: n, info
  
  ! 1 and 2
  integer, dimension(3) :: ipiv
  real(kind=8), dimension(3,3) :: a, alu, alu2, ain, work, c, q
  real(kind=8), dimension(3) :: b, x
  
  ! 3
  real(kind=8), dimension(8) :: work2
  real(kind=8), dimension(2,2) :: aa, vl, vr
  real(kind=8), dimension(2) :: wr, wi
  real(kind=8) :: alpha, beta
  
  ! 4
  real(kind=8), dimension(12) :: work3
  real(kind=8), dimension(3,3) :: vl3, vr3
  real(kind=8), dimension(3) :: wr3, wi3
  
  ! 5
  integer :: i, j
  real(kind=8), dimension(100,100) :: aaa, iiipiv
  real(kind=8), dimension(100) :: bbb
  
  print *, "==================="
  print *, "         1"
  print *, "==================="
  n = 3
  a = reshape([4.0d0,3.0d0,2.0d0,-2.0d0,6.0d0,1.0d0,1.0d0,-4.0d0,8.0d0],[3,3])
  alu = a
  call dgetrf(n,n,alu,n,ipiv,info)
  ain = alu
  call dgetri(n,ain,n,ipiv,work,n,info)
  call dgemm('N','N',n,n,n,1.0d0,a,n,ain,n,1.0d0,c,3)
  call dgemm('N','N',n,n,n,1.0d0,ain,n,a,n,1.0d0,q,3)
  
  print *, "Matrix A"
  print "(3f5.1)", transpose(a)
  print *, ""
  print *, "Inverse"
  print "(3f8.4)", transpose(ain)
  print *, ""
  print *, "Inverse*263"
  print "(3f5.1)", transpose(ain*263.0d0)
  print *, ""
  print *, "A*A^-1"
  print "(3f5.1)", transpose(c)
  print *, ""
  print *, "A^-1*A (DGEMM)"
  print "(3f5.1)", transpose(q)
  print *, "A^-1*A (MATMUL)"
  print "(3f5.1)", transpose(matmul(ain,a))
  print *, ""
  
  print *, "==================="
  print *, "         2"
  print *, "==================="
  b = (/12.0d0,-25.0d0,32.0d0/)
  alu = a
  x = b
  call dgesv(n,1,alu,n,ipiv,x,n,info)
  print "(a,3f7.1)", "b1 =", b
  print "(a,3f7.1)", "x1 =", x
  print *, ""

  b = (/4.0d0,-10.0d0,22.0d0/)
  alu = a
  x = b
  call dgesv(n,1,alu,n,ipiv,x,n,info)
  print "(a,3f7.1)", "b2 =", b
  print "(a,3f7.3)", "x2 =", x
  print *, ""

  b = (/20.0d0,-30.0d0,40.0d0/)
  alu = a
  x = b
  call dgesv(n,1,alu,n,ipiv,x,n,info)
  print "(a,3f7.1)", "b3 =", b
  print "(a,3f7.3)", "x3 =", x
  print *, ""
  
  print *, "==================="
  print *, "         3"
  print *, "==================="
  n = 2
  alpha = 1.0d0
  beta = 2.0d0
  aa = reshape([alpha,-beta,beta,alpha],[2,2])
  call dgeev('N','V',n,aa,n,wr,wi,vl,n,vr,n,work2,8,info)
  print *, "A"
  print "(2f5.1)", transpose(aa)
  print *, ""
  print *, "Real(lambda)"
  print "(2f7.3)", wr
  print *, "Imaginary(lambda)"
  print "(2f7.3)", wi
  print *, ""
  print *, "Eigenvectors"
  print *, " real, imaginary"
  print "(2f7.3)", transpose(vr)
  print *, ""
  
  print *, "==================="
  print *, "         4"
  print *, "==================="
  n = 3
  a = reshape([-2.0d0,2.0d0,-1.0d0,2.0d0,1.0d0,-2.0d0,-3.0d0,-6.0d0,0.0d0],[3,3])
  alu = a
  call dgeev('N','V',n,alu,n,wr3,wi3,vl3,n,vr3,n,work3,12,info)
  print *, "A"
  print "(3f5.1)", transpose(a)
  print *, ""
  print *, "Lambda"
  print "(3f5.1)", wr3
  print *, ""
  print *, "Eigenvectors"
  print *, "  x1     x2     x3"
  print "(3f7.3)", transpose(vr3)
  print *, ""
  
  print *, "==================="
  print *, "         5"
  print *, "==================="
  n = 100
  do i = 1,n
     do j = 1,n
        aaa(j,i) = 1.0d0/(i+j-1)
     enddo
  enddo
  do i = 1,n
     bbb(i)=1.0d0/i
  enddo
  call dgesv(n,1,aaa,n,iiipiv,bbb,n,info)
  print "(100f5.1)", bbb

end program main
