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

  integer,parameter :: L=501                  ! lattice side (square)
  integer,parameter :: d=2,DD=50
  integer,parameter :: n=(L-1)/(2*DD)
  integer,parameter,dimension(3) :: Lx=(/2,4,20/)
  real(kind=8),parameter :: eps=1.0d-6        ! convergence tolerance
  real(kind=8),parameter :: omega=1.80d0      ! SOR parameter
  real(kind=8),allocatable :: u(:,:)
  real(kind=8),allocatable :: chargedensity(:)
  real(kind=8) :: err,min,max,scale,charge
  integer :: i,j,k,m
  character(len=8) :: numstr
  
  allocate(u(L,L))
  
  do m=1,3
     ! initialize u to zero
     u=0.0d0
     
     ! make a plate on one side
     do j=L/2+1-Lx(m)*n,L/2+1+Lx(m)*n
        u(L/2+1+d/2*n,j)=0.5d0
     enddo
     
     ! Gauss-Seidel
     k=1
     do 
        !$omp parallel shared(u)
        !$omp do
        do j=2,L-1
           do i=L/2+2+mod(j,2),L-1,2
              if (.not.((i.eq.(L/2+1+n*d/2).and.((j.ge.(L/2+1-n*Lx(m)).and.(j.le.(L/2+1+n*Lx(m)))))))) then
                 u(i,j)=omega*0.25d0*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))+(1.0d0-omega)*u(i,j)
              endif
           enddo
        enddo
        !$omp do
        do j=2,L-1
           do i=L/2+2+mod(j+1,2),L-1,2
              if (.not.((i.eq.(L/2+1+n*d/2).and.((j.ge.(L/2+1-n*Lx(m)).and.(j.le.(L/2+1+n*Lx(m)))))))) then
                 u(i,j)=omega*0.25d0*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))+(1.0d0-omega)*u(i,j)
              endif
           enddo
        enddo
        !$omp end parallel
        
        ! calculate error estimate, every 100 iterations to
        ! save time
        if (mod(k,100)==0) then
           err=0.0d0
           !$omp parallel
           !$omp do reduction(+:err)
           do j=2,L-1
              do i=L/2+2,L-1
                 if (.not.((i.eq.(L/2+1+n*d/2).and.((j.ge.(L/2+1-n*Lx(m)).and.(j.le.(L/2+1+n*Lx(m)))))))) then
                    err=err-4.0d0*u(i,j)+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
                 endif
              enddo
           enddo
           !$omp end parallel
           if (abs(err)<eps) exit
        endif
        k=k+1
     enddo
     do i=L-1,L/2+2,-1
        u(L-i+1,:)=-u(i,:)
     enddo
     
     allocate(chargedensity(Lx(m)*2*n+1))
     write(numstr,'(i0.1)') m
     !top plate, upper surface
     open(21,file='topupper'//trim(numstr)//'.txt',status='replace')
     do j=1,2*Lx(m)*n+1
        chargedensity(j) = -real(n,8)*(u(L/2+1+d/2*n+1,L/2+1-Lx(m)*n+j-1)-u(L/2+1+d/2*n,L/2+1-Lx(m)*n+j-1))
        write(21,"(f10.6)") chargedensity(j)!-real(n,8)*(u(L/2+1+d/2*n+1,j)-u(L/2+1+d/2*n+1,j))
     enddo
     close(21)
     
     !top plate, lower surface
     open(21,file='toplower'//trim(numstr)//'.txt',status='replace')
     do j=1,2*Lx(m)*n+1
        chargedensity(j) = -real(n,8)*(u(L/2+1+d/2*n,L/2+1-Lx(m)*n+j-1)-u(L/2+1+d/2*n-1,L/2+1-Lx(m)*n+j-1))
        write(21,"(f10.6)") chargedensity(j)!-real(n,8)*(u(L/2+1+d/2*n+1,j)-u(L/2+1+d/2*n+1,j))
     enddo
     close(21)
     
     charge=0.0d0
     do j=1,2*Lx(m)*n
        charge = charge-(0.5d0)*(1.0d0/(real(n,8)))*(chargedensity(j)+chargedensity(j+1))
        !((u(L/2+1+d/2*n,j)-u(L/2+1+d/2*n-1,j))+((u(L/2+1+d/2*n,j+1)-u(L/2+1+d/2*n-1,j+1))))
     enddo
     write(*,*) 'Error ',abs(err),' in ',k,' iterations. c=C/L ',abs(charge/2), 'c/c0 ', abs(charge/real(Lx(m),8))
     deallocate(chargedensity)
     
     ! for testing speedup, skip writing out data below
     !stop
     
     ! write out solution
     open(20,file='laplace'//trim(numstr)//'.txt',status='replace')
     do i=1,L
        do j=1,L
           write(20,"(f10.6)",advance='no') u(i,j)
        enddo
        write(20,*)
     enddo
     close(20)
  enddo
  deallocate(u)

end program gseidel
