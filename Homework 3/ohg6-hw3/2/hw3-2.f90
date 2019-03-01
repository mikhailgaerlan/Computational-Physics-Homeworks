program differentiate
  implicit none
  
  integer :: i
  real(kind=8) :: d
  real(kind=8), dimension(3) :: x
  real(kind=8), dimension(3) :: answer
  
  x = (/0.1_8, 10.0_8, 100.0_8/)
  
  open(10, File='2/hw3-2-sin.txt', Status='Unknown')
  open(11, File='2/hw3-2-e.txt', Status='Unknown')
  
  do i = 1,3
     d = -2.0_8
     do while (d .le. 15.0_8)
        answer = fdiff(x(i),0.1_8*(10.0_8**(-d)))
        write(10, *) x(i), 0.1_8*(10.0_8**(-d)), answer, abs((answer-cos(x(i)))/(cos(x(i))))
        d=d+0.5
     end do
  end do
  
  answer=(/0.0_8,0.0_8,0.0_8/)
  
  do i = 1,3
     d = -2.0_8
     do while (d .le. 15.0_8)
        answer = gdiff(x(i),0.1_8*(10.0_8**(-d)))
        write(11, *) x(i), 0.1_8*(10.0_8**(-d)), answer, abs((answer-exp(x(i)))/(exp(x(i))))
        d=d+0.5
     end do
  end do
  
contains
  
  function f(x)
    real(kind=8) :: f, x
    
    f = sin(x)
    
  end function f

  function g(x)
    real(kind=8) :: g, x
    
    g = exp(x)
    
  end function g

  function fdiff(xn, h)
    real(kind=8) :: xn, h
    real(kind=8), dimension(-2:2) :: fv
    real(kind=8), dimension(3) :: fdiff
    integer :: i
    
    do i=-2,2
       fv(i)=f(xn+i*h)
    end do
    fdiff(1) = (fv(2)-fv(0))/(2*h)
    fdiff(2) = (fv(1)-fv(-1))/(2*h)
    fdiff(3) = (fv(-2)-8*fv(-1)+8*fv(1)-fv(2))/(12*h)
    
  end function fdiff

  function gdiff(xn, h)
    real(kind=8) :: xn, h
    real(kind=8), dimension(-2:2) :: fv
    real(kind=8), dimension(3) :: gdiff
    integer :: i
    
    do i=-2,2
       fv(i)=g(xn+i*h)
    end do
    gdiff(1) = (fv(2)-fv(0))/(2*h)
    gdiff(2) = (fv(1)-fv(-1))/(2*h)
    gdiff(3) = (fv(-2)-8*fv(-1)+8*fv(1)-fv(2))/(12*h)
    
  end function gdiff
  
end program differentiate
