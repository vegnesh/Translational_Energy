module BINVEL
CONTAINS
subroutine bin_velocity(velocity_value,beta,c_initial,c_final,w)
 USE global_parameters, only  : PI
 USE basic_functions
 IMPLICIT NONE
 REAL,DIMENSION(3),intent(in) ::c_initial,c_final
 REAL,DIMENSION(3),intent(out) :: velocity_value
 REAL,DIMENSION(3) :: wbin
 REAL,DIMENSION(3),optional :: w
 REAL,intent(in):: beta
 INTEGER :: i
  if(present(w)) then 
 wbin = w
 else
 wbin(1) = 0.D00
 wbin(2) = 0.D00
 wbin(3) = 0.D00
end if 
do i=1,3 
 velocity_value(i) =wbin(i) + (1/sqrt(PI*beta))&
		   *(binexp(beta,c_initial(i),wbin(i))&
		   - binexp(beta,c_final(i),wbin(i)))/ &
		   (binerf(beta,c_final(i),wbin(i))&
		   -binerf(beta,c_initial(i),wbin(i)))
end do
end subroutine bin_velocity 
end module BINVEL
!***********************************************************!

module BINENERGY
CONTAINS
subroutine bin_energy(energy_value,beta,c_initial,c_final,w)
USE global_parameters, only : PI
USE BINVEL
USE basic_functions
 IMPLICIT NONE
 REAL,DIMENSION(3),intent(in) ::c_initial,c_final
 REAL,intent(out) :: energy_value
 REAL,DIMENSION(3) :: wbin
 REAL,DIMENSION(3),optional :: w
 REAL,intent(in):: beta
 INTEGER :: i
  REAL,DIMENSION(3) :: tempvar
 
 if(present(w)) then 
 wbin = w
 else
 wbin(1) = 0.D00
 wbin(2) = 0.D00
 wbin(3) = 0.D00
end if
  
 call bin_velocity(tempvar,beta,c_initial,c_final,wbin)
! print *, "temp value" , tempvar
 energy_value = DOT_PRODUCT(wbin,wbin) + DOT_PRODUCT(wbin,tempvar-wbin) + 1.5/beta + &
                    (1/sqrt(PI*beta))&
		   *(c_initial(1)*binexp(beta,c_initial(1),wbin(1))&
		   - c_final(1)*binexp(beta,c_final(1),wbin(1)))/ &
		   (binerf(beta,c_final(1),wbin(1))&
		   -binerf(beta,c_initial(1),wbin(1)))&
	         +(1/sqrt(PI*beta))&
		   *(c_initial(2)*binexp(beta,c_initial(2),wbin(2))&
		   - c_final(2)*binexp(beta,c_final(2),wbin(2)))/ &
		   (binerf(beta,c_final(2),wbin(2))&
		   -binerf(beta,c_initial(2),wbin(2)))&
	         +(1/sqrt(PI*beta))&
		   *(c_initial(3)*binexp(beta,c_initial(3),wbin(3))&
		   - c_final(3)*binexp(beta,c_final(3),wbin(3)))/ &
		   (binerf(beta,c_final(3),wbin(3))&
		   -binerf(beta,c_initial(3),wbin(3)))
 
 
end subroutine bin_energy
end module BINENERGY
