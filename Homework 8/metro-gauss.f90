!
! Example: use Metropolis method to generate a Gaussian distribution
! exp(-(x-x_0)**2) in one dimension
! 
! RT Clay  10/2015
!
program metro1d
  implicit none

  integer,parameter :: neq=500             ! equilibration steps
  integer,parameter :: nmeas=1000000       ! measurements steps
  real(kind=8),parameter :: stepsize=1.0d0 ! Metropolis step size (delta_x)

  ! for making a binned plot of the results
  real(kind=8),parameter :: xmax=10.0d0   ! maximum x on p(x) plot
  integer,parameter :: nbin=100           ! number of bins in interval (0,xmax)
  integer,dimension(nbin) :: bin

  real(kind=8) :: x,temp,weight
  integer :: i,j,nacc

  ! initialize random generator
  call set_seed
        
  ! choose a starting configuration: random x between 0 and xmax
  call random_number(x)
  x=x*xmax
  weight=w(x)

  ! file to print out series of x's
  open(unit=12,file='mc.dat',status='replace')

  ! equilibrate
  do i=1,neq
     call mc_step(x,weight,stepsize,nacc)
     write(12,*) i,x
  enddo
  
  ! measure
  bin=0
  nacc=0

  do i=1,nmeas
     call mc_step(x,weight,stepsize,nacc)
     write(12,*) i+neq,x

     ! x should have Gaussian distribution. Make a histogram
     ! of the results
     temp=x*nbin/xmax
     j=temp+1
     if (j<nbin) bin(j)=bin(j)+1
  enddo
  close(12)

  ! normalize and write out p(x)
  open(unit=12,file='gauss.dat',status='replace')
  do i=1,nbin
     write(12,*) (xmax/nbin*0.5d0)+(i-1)*xmax/nbin,real(bin(i))/nmeas*nbin/xmax
  enddo
  close(12)

  ! compute acceptance rate
  write(*,*) 'Metropolis finished. Acceptance rate=',nacc/real(nmeas)


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
    real(kind=8) :: x,r
    integer,intent(inout) :: nacc

    logical :: iacc
    real(kind=8) :: deltax,neww,xx,temp
    real(kind=8) :: newx
      
    ! make random change in x between [-stepsize,stepsize]
    call random_number(xx)
    deltax=stepsize*2.0d0*(0.5d0-xx)
    temp=deltax+x

    ! check to be sure x stays within limits of integration (0,xmax)
    if (temp<0.0d0.or.temp>xmax) then
       newx=x
       ! don't have to do anything else, exit subroutine
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
       call random_number(temp)
       if (r>=temp) then
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
    real(kind=8),intent(in) :: x
    real(kind=8) :: w
    real(kind=8) :: sum
    integer :: i

    w=exp(-(x-xmax*0.5d0)**2.0d0)

  end function w

  !
  ! set fortran random seed from system clock
  !
  subroutine set_seed
    integer,allocatable :: seed(:)
    integer :: n_seed,clock,i
    real(kind=8) :: x
    
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


end program metro1d
