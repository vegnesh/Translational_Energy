program test
use jacobian
use velocityjacobian
use global_parameters, only :PI
use binenergy
use binvel
use basic_functions
use matinv
IMPLICIT NONE
REAL(KIND=8) :: x,y,z,h,betaguess,eval
external RROR
REAL,DIMENSION(3) :: c_initial,c_final,wbin,fgh,wbin2,velval
REAL, DIMENSION(4,4) :: A,B,C
REAL, DIMENSION(4,1) :: rhsvec,xnew,xold
INTEGER, PARAMETER :: seed = 86456
INTEGER :: i,j,k
fgh(1) = 5.0
fgh(2) = 6.0
fgh(3) = 8.0

write(*, '(A)', ADVANCE = "NO") "Enter the value of beta:  "
read(*,*) x
write(*, '(A)', ADVANCE = "NO") "Enter the value of Min Velocity:  "
read(*,*) c_initial
write(*, '(A)', ADVANCE = "NO") "Enter the value of Max Velocity:  "
read(*,*) c_final
write(*, '(A)', ADVANCE = "NO") "Enter the value of Bin Velocity:  "
read(*,*) wbin


h = exp(x)
call error(x,z)
wbin2 = wbin
wbin2(1) = wbin2(1)*1.0000001
call bin_velocity(fgh,x,c_initial,c_final,wbin)
call bin_energy(h,x,c_initial,c_final,wbin)
call bin_energy(y,x,c_initial,c_final,wbin2)



print *,wbin

print *,"The value of BINVEL is :",fgh
print *,"The value of BINENERGY is :",h
print *,"The value of the Energy Jacobian for beta is :",jacenergy_beta(x,c_initial,c_final,wbin)
print *,"The value of Energy Jacobian along x is:",jacenergy_w(x,c_initial(1),c_final(1),wbin(1)),(y-h)/(wbin2(1)-wbin(1))
print *,"The value of Energy Jacobian along y is:",jacenergy_w(x,c_initial(2),c_final(2),wbin(2))
print *,"The value of Energy Jacobian along z is:",jacenergy_w(x,c_initial(3),c_final(3),wbin(3))
print *,"The value of the Velocity Jacobian along x for beta is :",jacvel_beta(x,c_initial(1),c_final(1),wbin(1))
print *,"The value of the Velocity Jacobian along y for beta is :",jacvel_beta(x,c_initial(2),c_final(2),wbin(2))
print *,"The value of the Velocity Jacobian along z for beta is :",jacvel_beta(x,c_initial(3),c_final(3),wbin(3))
print *,"The value of Velocity Jacobian along x is:",jacvel_w(x,c_initial(1),c_final(1),wbin(1))
print *,"The value of Velocity Jacobian along y is:",jacvel_w(x,c_initial(2),c_final(2),wbin(2))
print *,"The value of Velocity Jacobian along z is:",jacvel_w(x,c_initial(3),c_final(3),wbin(3))

do i=1,4
do j=1,4
 
 A(i,j) = 0.D0
 
end do
end do
A(1,1) = jacenergy_beta(x,c_initial,c_final,wbin)
A(2,1) = jacvel_beta(x,c_initial(1),c_final(1),wbin(1))
A(3,1) = jacvel_beta(x,c_initial(2),c_final(2),wbin(2))
A(4,1) = jacvel_beta(x,c_initial(3),c_final(3),wbin(3))
A(1,2) = jacenergy_w(x,c_initial(1),c_final(1),wbin(1))
A(1,3) = jacenergy_w(x,c_initial(2),c_final(2),wbin(2)) 
A(1,4) = jacenergy_w(x,c_initial(3),c_final(3),wbin(3))
A(2,2) = jacvel_w(x,c_initial(1),c_final(1),wbin(1)) 
A(3,3) = jacvel_w(x,c_initial(2),c_final(2),wbin(2))
A(4,4) = jacvel_w(x,c_initial(3),c_final(3),wbin(3))

B = matinv4(A)

write(*, '(A)', ADVANCE = "NO") "Enter the value of beta guess:  "
read(*,*) betaguess
write(*, '(A)', ADVANCE = "NO") "Enter the value of Wguess:  "
read(*,*) wbin2

do k=1,5
   call bin_velocity(velval,betaguess,c_initial,c_final,wbin2)
   call bin_energy(eval,betaguess,c_initial,c_final,wbin2)
   do i=1,4
   do j=1,4
    
    A(i,j) = 0.D0
    
   end do
   end do
   A(1,1) = jacenergy_beta(betaguess,c_initial,c_final,wbin2)
   A(2,1) = jacvel_beta(betaguess,c_initial(1),c_final(1),wbin2(1))
   A(3,1) = jacvel_beta(betaguess,c_initial(2),c_final(2),wbin2(2))
   A(4,1) = jacvel_beta(betaguess,c_initial(3),c_final(3),wbin2(3))
   A(1,2) = jacenergy_w(betaguess,c_initial(1),c_final(1),wbin2(1))
   A(1,3) = jacenergy_w(betaguess,c_initial(2),c_final(2),wbin2(2)) 
   A(1,4) = jacenergy_w(betaguess,c_initial(3),c_final(3),wbin2(3))
   A(2,2) = jacvel_w(betaguess,c_initial(1),c_final(1),wbin2(1)) 
   A(3,3) = jacvel_w(betaguess,c_initial(2),c_final(2),wbin2(2))
   A(4,4) = jacvel_w(betaguess,c_initial(3),c_final(3),wbin2(3))
   B = matinv4(A)
   
   rhsvec(1,1) = eval
   rhsvec(2:4,1) = velval
   xold(1,1) = betaguess
   xold(2:4,1) = wbin2
   
   xnew = xold - matmul(B,rhsvec)
   betaguess = xnew(1,1)
   wbin2 = xnew(2:4,1)
print *,betaguess , wbin2

end do

!
!print *, "Matrix A"
!do i=1,4
!  print *, A(i,1),A(i,2),A(i,3),A(i,4)
!end do
!print *, "Matrix B"
!do i=1,4
!  print *, B(i,1),B(i,2),B(i,3),B(i,4)
!end do
!print *, "Matrix C"
!C = matmul(A,B)
!do i=1,4
!  print *, C(i,1),C(i,2),C(i,3),C(i,4)
!end do



end program test
