program homework
  implicit none

  integer :: i, n
  real(kind=4) :: incsum, decsum
  real(kind=4), parameter :: ans = 1.2337005501361698273543113749845188919142124259050988
  
  print *, "       nmax     inc. sum         dec. sum           inc. error       dec. error    conclusion"
  print *, "---------------------------------------------------------------------------------------------"
  do n=2000,8000,1000
     incsum=0.0
     decsum=0.0
     do i=0,n
        incsum = incsum + 1.0/real((2*i+1)**2)
        decsum = decsum + 1.0/real((2*(n-i)+1)**2)
     enddo
     if (incsum==decsum) then
        print *, n, incsum, decsum, abs(incsum-ans)/ans, abs(decsum-ans)/ans, "   EQUAL"
     else if (incsum-ans < decsum-ans) then
        print *, n, incsum, decsum, abs(incsum-ans)/ans, abs(decsum-ans)/ans, "   DECREASING"
     else
        print *, n, incsum, decsum, abs(incsum-ans)/ans, abs(decsum-ans)/ans, "   INCREASING"
     endif
  enddo
end program homework
