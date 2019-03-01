program homework
  implicit none

  integer :: i, n, io, iu
  integer :: xi, underi, overi
  real(kind=4) :: x, under, over
  real(kind=8) :: x2, under2, over2

  print *, "Single Precision"
  print *, " n     over   under"
  print *, "-------------------"  
  do n=10,1,-1
     under = 1.0
     over = 1.0
     
     i=0
     io=0
     iu=0
     x = real(n/10.0+1.0)
     do
        i=i+1
        under = under / x
        over = over * x
        
        !write (*,*) i, over, under
        if (over > huge(over) .and. io==0) then
           io=i
        else if (under==under/x .and. iu==0) then
           iu=i
        else if ((io .ne. 0) .and. (iu .ne. 0)) then
           exit
        end if
     end do
     print "(a,f3.1,a,i4,a,i4)", " ", x,"    ", io,"    ", iu
  end do
  
  print *, ""
  print *, "Double Precision"
  print *, " n     over   under"
  print *, "-------------------"
  do n=10,1,-1
     under2 = 1.0_8
     over2 = 1.0_8
     
     i=0
     io=0
     iu=0
     x2 = real(n/10.0_8+1.0_8,8)
     do
        i=i+1
        under2 = under2 / x2
        over2 = over2 * x2
        
        !write (*,*) i, over, under
        if (over2 > huge(over2) .and. io==0) then
           io=i
        else if (under2==under2/x2 .and. iu==0) then
           iu=i
        else if ((io .ne. 0) .and. (iu .ne. 0)) then
           exit
        end if
     end do
     print "(a,f3.1,a,i4,a,i4)", " ", x2,"    ", io,"    ", iu
  end do

  print *, ""
  print *, "Integer"
  print *, "-------------------"
  underi = 0
  overi = 0
  
  i=0
  io=0
  iu=0
  do
     underi = underi - 1
     overi = overi + 1
     if(real(overi)/real(overi+1)<0)exit
  end do
  print *, overi, underi

end program homework
