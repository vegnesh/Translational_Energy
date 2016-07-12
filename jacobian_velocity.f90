module velocityjacobian
CONTAINS
Real function jacvel_w(beta,c_initial,c_final,w) 
   USE global_parameters, only : PI,K_B
   USE basic_functions
   IMPLICIT NONE
   REAL,INTENT(in) :: beta
   REAL,INTENT(in) :: w
   REAL,INTENT(in) :: c_initial,c_final
   REAL :: NUMER, DENO, DNUMER,DDENO
   !REAL :: beta
   INTEGER :: i = 0
   jacvel_w = 0.D0
   ! Term 1
   jacvel_w = jacvel_w + 1.D0

   ! Term 2 
   NUMER = (binexp(beta,c_initial,w)&
   	   -binexp(beta,c_final,w))

   DENO  =  (binerf(beta,c_final,w)&
            - binerf(beta,c_initial,w))
  
   DDENO = -2.0*sqrt(beta/PI) * (binexp(beta,c_final,w)&
           - binexp(beta,c_initial,w))
   DNUMER = 2*beta*((c_initial-w)*binexp(beta,c_initial,w)&
            - (c_final-w)*binexp(beta,c_final,w))
   
   jacvel_w = jacvel_w +  (1/sqrt(PI*beta))&
                *(DNUMER*DENO - DDENO*NUMER)/DENO**2

end function jacvel_w

Real function jacvel_beta(beta,c_initial,c_final,w)
   USE global_parameters, only : PI,K_B
   USE basic_functions
   IMPLICIT NONE
   REAL,INTENT(in) :: beta
   REAL,INTENT(in) :: w
   REAL,INTENT(in) :: c_initial,c_final
   REAL :: NUMER, DENO, DNUMER,DDENO
   !REAL :: beta
   INTEGER :: i = 0
   jacvel_beta = 0.D0
     NUMER  = (binexp(beta,c_initial ,w )&
   	   -binexp(beta,c_final ,w ))

     DENO   = (binerf(beta,c_final ,w ) -&
                binerf(beta,c_initial ,w ))

     DDENO  = (1/(sqrt(PI*beta))) *&
               ((c_final -w )*binexp(beta,c_final ,w )&
   	- (c_initial -w )*binexp(beta,c_initial ,w ))

     DNUMER = ((-(c_initial -w )**2)*&
                binexp(beta,c_initial ,w )&
   	         +((c_final -w )**2)*&
                binexp(beta,c_final ,w ))
     jacvel_beta = jacvel_beta - (1/(sqrt(PI*beta)*2*beta))*&
             (NUMER /DENO )
     jacvel_beta = jacvel_beta +  (1/(sqrt(PI*beta)))*(DENO *DNUMER &
                 - NUMER *DDENO )/(DENO **2)


end function jacvel_beta

end module velocityjacobian

