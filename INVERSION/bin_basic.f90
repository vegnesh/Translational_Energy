!In this file the basic function required for defining the bin !
!parameters are defined!
module basic_functions
CONTAINS
REAL FUNCTION binexp(beta,c,w)
!Takes in values of bin parameters beta and w, the velocity c and returns!
! the value of the bin exponent !
 IMPLICIT NONE
 REAL :: beta,c,wval
 REAL,OPTIONAL :: w

 if(present(w)) then
  wval = w
 else 
  wval = 0.D00
 endif

 binexp = exp(-beta*((c-wval)**2))

END FUNCTION binexp

REAL FUNCTION binerf(beta,c,w)
!Takes in values of bin parameters beta and w, the velocity c and returns!
! the value of the bin error !

!NOTE BETA SHOULD NOT BE NEGATIVE!

 IMPLICIT NONE
 REAL ::  beta,c,wval,temp1,temp2
 REAL,OPTIONAL :: w
 EXTERNAL :: ERROR
 if(present(w)) then
  wval = w
 else
  wval = 0.D00
 endif
 
 temp1 = (c-wval)*sqrt(beta)
 call ERROR(temp1,temp2)
 binerf = temp2

END FUNCTION binerf
end module basic_functions
