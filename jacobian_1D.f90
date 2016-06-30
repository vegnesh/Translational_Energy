module jacobian
CONTAINS
Real FUNCTION jac1D(beta,c_initial,c_final,w)
   USE global_parameters, only : PI,K_B
   USE basic_functions
   IMPLICIT NONE
   REAL,INTENT(in) :: beta
   REAL,INTENT(in),optional,DIMENSION(3) :: w
   REAL,DIMENSION(3),INTENT(in) :: c_initial,c_final
   REAL,DIMENSION(3) :: wbin
      REAL, DIMENSION(3) :: NUMER, DENO, DNUMER,DDENO
   !REAL :: beta
   INTEGER :: i = 0
   if (present(w)) then 
    wbin = w
   else 
    wbin(1) = 0.D0
    wbin(2) = 0.D0
    wbin(3) = 0.D0
   end if
   print *,w, wbin, "c:" , c_initial
   jac1D = 0.D00
   !Term 4
   jac1D = jac1D - 1.5/(beta)**2 
   
   !Term 2
   do i=1,3
     NUMER(i) = wbin(i)*(binexp(beta,c_initial(i),wbin(i))&
   	   -binexp(beta,c_final(i),wbin(i)))

     DENO(i)  = (binerf(beta,c_final(i),wbin(i)) -&
                binerf(beta,c_initial(i),wbin(i)))

     DDENO(i) = (1/(sqrt(PI*beta))) *&
               ((c_final(i)-wbin(i))*binexp(beta,c_final(i),wbin(i))&
   	- (c_initial(i)-wbin(i))*binexp(beta,c_initial(i),wbin(i)))

     DNUMER(i)= wbin(i)*((-(c_initial(i)-wbin(i))**2)*&
                binexp(beta,c_initial(i),wbin(i))&
   	         +((c_final(i)-wbin(i))**2)*&
                binexp(beta,c_final(i),wbin(i)))
     jac1D = jac1D - (1/(sqrt(PI*beta)*2*beta))*&
             (NUMER(i)/DENO(i))
     jac1D = jac1D +  (1/(sqrt(PI*beta)))*(DENO(i)*DNUMER(i)&
                 - NUMER(i)*DDENO(i))/(DENO(i)**2)


   end do
   
   !Term 3 first part
   do i=1,3
     jac1D = jac1D -  (1/(sqrt(PI*beta)*2*beta))*&
     (c_initial(i)*binexp(beta,c_initial(i),wbin(i))&
     - c_final(i)*binexp(beta,c_final(i),wbin(i)))/ &
     (binerf(beta,c_final(i),wbin(i))&
     -binerf(beta,c_initial(i),wbin(i)))
   end do
   !Term 3 Part 2
   do i=1,3
    NUMER(i) = (c_initial(i)*binexp(beta,c_initial(i),wbin(i))&
   	   -c_final(i)*binexp(beta,c_final(i),wbin(i)))
    
    DENO(i) = (binerf(beta,c_final(i),wbin(i)) -&
              binerf(beta,c_initial(i),wbin(i)))
   
    DNUMER(i) =  (-c_initial(i)*(c_initial(i)-wbin(i))**2)*&
                binexp(beta,c_initial(i),wbin(i))&
   	         +(c_final(i)*(c_final(i)-wbin(i))**2)*&
                binexp(beta,c_final(i),wbin(i))
    DDENO(i)  = (1/(sqrt(PI*beta))) *&
        ((c_final(i)-wbin(i))*binexp(beta,c_final(i),wbin(i))&
   	- (c_initial(i)-wbin(i))*binexp(beta,c_initial(i),wbin(i)))
               
    jac1D = jac1D +  (1/(sqrt(PI*beta)))*(DENO(i)*DNUMER(i)&
                 - NUMER(i)*DDENO(i))/(DENO(i)**2)
  !  jac1D = jac1D +  (1/(sqrt(PI*beta)))*(DENO(i)&
  ! 	*(-(c_initial(i)**3)*exp(-beta*c_initial(i)**2)&
  ! 	+(c_final(i)**3)*exp(-beta*c_final(i)**2))&
  !      - NUMER(i)*(1/(sqrt(PI*beta)))*(c_final(i)*exp(-beta*c_final(i)**2)&
  ! 	- c_initial(i)*exp(-beta*c_initial(i)**2)))/(DENO(i)**2)

   end do
  
end FUNCTION jac1D
end module jacobian

!   jac1D = jac1D - (1/(sqrt(PI*beta)*2*beta))*&
!   		 ((c_initial(1)*binexp(beta,c_initial(1),wbin(1))&
!   		   - c_final(1)*binexp(beta,c_final(1),wbin(1)))/ &
!   		   (binerf(beta,c_final(1),wbin(1))&
!   		   -binerf(beta,c_initial(1),wbin(1)))&
!   	         +(c_initial(2)*binexp(beta,c_initial(2),wbin(2))&
!   		   - c_final(2)*binexp(beta,c_final(2),wbin(2)))/ &
!   		   (binerf(beta,c_final(2),wbin(2))&
!   		   -binerf(beta,c_initial(2),wbin(2)))&
!   	         +(c_initial(3)*binexp(beta,c_initial(3),wbin(3))&
!   		   - c_final(3)*binexp(beta,c_final(3),wbin(3)))/ &
!   		   (binerf(beta,c_final(3),wbin(3))&
!   		   -binerf(beta,c_initial(3),wbin(3))))
!   do i =1,3
!   NUMER(i) = (c_initial(i)*binexp(beta,c_initial(i),wbin(i))&
!   	   -c_final(i)*binexp(beta,c_final(i),wbin(i)))
!    
!   DENO(i) = (binerf(beta,c_final(i),wbin(i)) - binerf(beta,c_initial(i),wbin(i)))
!   
!   
!   
!   jac1D = jac1D + (1/(sqrt(PI*beta)))*(DENO(i)&
!   	*(-(c_initial(i)**3)*exp(-beta*c_initial(i)**2)&
!   	+(c_final(i)**3)*exp(-beta*c_final(i)**2))&
!   	 - NUMER(i)*(1/(sqrt(PI*beta)))*(c_final(i)*exp(-beta*c_final(i)**2)&
!   	- c_initial(i)*exp(-beta*c_initial(i)**2)))/(DENO(i)**2)
!   end do


!******************************************************
!    DNUMER(i) =  (-c_initial(i)*(c_initial(i)-wbin(i))**2)*&
!                binexp(beta,c_initial(i),wbin(i))&
!   	         +(c_final(i)*(c_final(i)-wbin(i))**2)*&
!                binexp(beta,c_final(i),wbin(i))
!    DDENO(i)  = (1/(sqrt(PI*beta))) *&
!        ((c_final(i)-wbin(i))*binexp(beta,c_final(i),wbin(i))&
!   	- (c_initial(i)-wbin(i))*binexp(beta,c_initial(i),wbin(i)))
                
    !jac1D = jac1D +  (1/(sqrt(PI*beta)))*(DENO(i)*DNUMER(i)&
    !              - NUMER(i)*DDENO(i))/(DENO(i)**2)
!    jac1D = jac1D +  (1/(sqrt(PI*beta)))*(DENO(i)&
!   	*(-(c_initial(i)**3)*exp(-beta*c_initial(i)**2)&
!   	+(c_final(i)**3)*exp(-beta*c_final(i)**2))&
!        - NUMER(i)*(1/(sqrt(PI*beta)))*(c_final(i)*exp(-beta*c_final(i)**2)&
!   	- c_initial(i)*exp(-beta*c_initial(i)**2)))/(DENO(i)**2)

! 
