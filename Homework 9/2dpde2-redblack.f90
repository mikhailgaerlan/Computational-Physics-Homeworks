! 2D PDE solution via Gauss-Seidel relaxation
!
!  gfortran -O2 -o 2dpde.exe 2dpde2-redblack.f90
!
! uses red/black decomposition
!
! RT Clay 11/2015
!
program gseidel
  implicit none

  integer,parameter :: L=256                  ! lattice side (square)
  real(kind=8),parameter :: eps=1.0d-6        ! convergence tolerance
  real(kind=8),parameter :: omega=1.80d0      ! SOR parameter
  real(kind=8),allocatable :: u(:,:)
  real(kind=8) :: err,min,max,scale
  integer :: i,j,k,n

  allocate(u(L,L))

  ! initialize u to zero
  u=0.0d0

  ! make a plate on one side
  do i=L/2-50,L/2+50
    u(i,1)=0.5d0
  enddo
  
  ! Gauss-Seidel
  k=1
  do 
     do j=2,L-1
        do i=2+mod(j,2),L-1,2
           u(i,j)=omega*0.25d0*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))+(1.0d0-omega)*u(i,j)
        enddo
     enddo
     do j=2,L-1
        do i=2+mod(j+1,2),L-1,2
           u(i,j)=omega*0.25d0*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))+(1.0d0-omega)*u(i,j)
        enddo
     enddo

     ! calculate error estimate, every 100 iterations to
     ! save time
     if (mod(k,100)==0) then
        err=0.0d0
        do j=2,L-1
           do i=2,L-1
              err=err-4.0d0*u(i,j)+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
           enddo
        enddo
        
        if (abs(err)<eps) exit
     endif
     k=k+1
  enddo
  write(*,*) 'Error ',abs(err),' in ',k,' iterations.'

  ! for testing speedup, skip writing out data below
  !stop

  ! write out solution
  open(20,file='laplace.dat',status='replace')
  do i=1,L
     do j=1,L
        write(20,"(f10.6)",advance='no') u(i,j)
     enddo
     write(20,*)
  enddo
  close(20)
  deallocate(u)

end program gseidel
