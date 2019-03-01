program homework
  implicit none
  
  integer:: i
  real(kind=4) :: eps1, one1
  real(kind=8) :: eps2, one2
  
  print *, "Single Precision"
  print *, "----------------"
  eps1 = 1.0
  one1 = 1.0
  do
     eps1 = eps1 / 2.0
     one1 = 1.0 + eps1
     if (1.0 == one1) then
        write (*,*) eps1
        exit
     endif
  end do

  print *, ""
  print *, "Double Precision"
  print *, "----------------"
  eps2 = 1.0_8
  one2 = 1.0_8
  do
     eps2 = eps2 / 2.0_8
     one2 = 1.0_8 + eps2
     if (1.0_8 == one2) then
        write (*,*) eps2
        exit
     endif
  end do

end program homework
