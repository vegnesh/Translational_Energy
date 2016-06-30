program ERFTEST
  IMPLICIT NONE
  INTEGER :: i
  INTEGER,parameter :: N = 7001
  REAL, DIMENSION(N) :: x, erfval
  EXTERNAL :: ERROR
  
  do i=1,N
     x(i) = -3.5 + (i-1)*0.001
     call ERROR(x(i),erfval(i))
     print *, x(i), erfval(i)
  end do

end program ERFTEST
