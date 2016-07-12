program test
use jacobian
use velocityjacobian
use global_parameters, only :PI
use binenergy
use binvel
use basic_functions
IMPLICIT NONE
REAL(KIND=8) :: x,y,z,h
external ERROR
REAL,DIMENSION(3) :: c_initial,c_final,wbin,fgh,wbin2
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



end program test
