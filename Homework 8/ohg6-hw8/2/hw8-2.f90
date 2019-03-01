!
! Example: use Metropolis method to generate a Gaussian distribution
! exp(-(x-x_0)**2) in one dimension
! 
! RT Clay  10/2015
!
program metro1d
  implicit none

  integer,parameter :: neq=500                       ! equilibration steps
  integer,parameter :: nmeas=10000000                 ! measurements steps
  real(kind=8),parameter :: stepsize=1.0d0           ! Metropolis step size (delta_x)

  ! for making a binned plot of the results
  real(kind=8),parameter :: xmax=5.0d0,xmin=-5.0d0   ! maximum x on p(x) plot
  integer,parameter :: nbin=100,dim=10               ! number of bins in interval (0,xmax)
  integer,dimension(nbin) :: bin
  
  real(kind=8) :: weight,temp
  real(kind=8),dimension(dim) :: x  
  integer :: i,j,nacc
  real(kind=8) :: tot, tot2, tot4
  
  ! initialize random generator
  call set_seed
        
  ! choose a starting configuration: random x between 0 and xmax
  call random_array(x)
  x=x*(xmax-xmin)+xmin
  weight=w(x)

  ! file to print out series of x's
  open(unit=12,file='mc.txt',status='replace')
  
  ! equilibrate
  do i=1,neq
     call mc_step(x,weight,stepsize,nacc)
     write(12,*) i,x
  enddo
  
  ! measure
  bin=0
  nacc=0
  tot=0.0d0
  tot2=0.0d0
  tot4=0.0d0
  
  do i=1,nmeas
     call mc_step(x,weight,stepsize,nacc)
     write(12,*) i+neq,x
     
     ! x should have Gaussian distribution. Make a histogram
     ! of the results
     tot=tot+x(1)
     tot2=tot2+x(1)**2
     tot4=tot4+x(1)**4
     temp=(x(1)-xmin)*nbin/(xmax-xmin)
     j=temp+1
     if (j<nbin) bin(j)=bin(j)+1
  enddo
  close(12)
  
  ! normalize and write out p(x)
  open(unit=12,file='gauss.txt',status='replace')
  do i=1,nbin
     write(12,*) xmin+((xmax-xmin)/nbin*0.5d0)+(i-1)*(xmax-xmin)/nbin,real(bin(i))/nmeas*nbin/(xmax-xmin)
  enddo
  close(12)
  
  ! compute acceptance rate
  write(*,*) 'Metropolis finished. Acceptance rate=',nacc/real(nmeas)
  write(*,*) 'Expectation value of <x>=',tot/nmeas
  write(*,*) 'Expectation value of <x^2>=',tot2/nmeas
  write(*,*) 'Expectation value of <x^4>=',tot4/nmeas

contains

  subroutine mc_step(x,weight,stepsize,nacc)
    implicit none
    ! take a single Metropolis step 
    !     
    !  x : position x
    !  weight : current weight
    !  stepsize : size of step
    !  nacc : counts number of accepted moves
    !
    real(kind=8),intent(in) :: stepsize
    real(kind=8),intent(inout) :: weight
    real(kind=8), dimension(dim) :: x
    integer,intent(inout) :: nacc
    integer :: i
    
    logical :: iacc,limits
    real(kind=8), dimension(dim) :: deltax,newx,xx,temp
    real(kind=8) :: r,neww,temp2
    
    ! make random change in x between [-stepsize,stepsize]
    call random_array(xx)
    deltax=stepsize*2.0d0*(0.5d0-xx)
    temp=deltax+x

    ! check to be sure x stays within limits of integration (xmin,xmax)
    limits=.false.
    do i=1,dim
       if (temp(i)<xmin.or.temp(i)>xmax) then
          limits=.true.
          exit
       endif
    enddo
    if (limits) then
       newx=x
       nacc=nacc+1
       return
    else
       newx=temp
    endif
    
    ! new weight
    neww=w(newx)

    ! Metropolis algorithm is here:
    ! compare ratio of weights to uniform random #
    r=neww/weight
    iacc=.false.

    ! check first to see if the ration is bigger than 1
    ! in this case, we do not need to generate a random #
    if (r>1.0d0) iacc=.true.

    ! ratio less than 1, need to compare to random #
    if (iacc.eqv..false.) then
       call random_number(temp2)
       if (r>=temp2) then
          iacc=.true.
       endif
    endif

    !  accept the move 
    if (iacc) then
       x=newx
       weight=neww
       nacc=nacc+1
    endif
  end subroutine mc_step

!
! weight function (not normalized) 
!
! Note: for this weight function, can compute ratio with only a 
! single call to exp(x)
!
  function w(x) 
    real(kind=8),intent(in),dimension(:) :: x
    real(kind=8) :: w

    w=product(exp(-0.5d0*(x**2.0d0)))

  end function w

  !
  ! set fortran random seed from system clock
  !
  subroutine set_seed
    integer,allocatable :: seed(:)
    integer :: n_seed,clock,i
    
    ! get the number of seed integers (varies with system)
    !
    ! note use of 'keyword argument' size
    call random_seed(size=n_seed)
    
    ! allocate seed array
    allocate(seed(n_seed))
    
    ! initialize the seed from the time
    call system_clock(count=clock)

    do i=1,n_seed
       seed(i)=clock/i  ! could do other things here
    enddo

    ! actually initialize the generator
    call random_seed(put=seed)

    deallocate(seed)
  end subroutine set_seed

  subroutine random_array(x)
    real(kind=8) :: rand
    real(kind=8), dimension(:) :: x
    integer :: i
    
    do i=1,dim
       call random_number(rand)
       x(i) = rand
    enddo
  end subroutine random_array
  
end program metro1d
