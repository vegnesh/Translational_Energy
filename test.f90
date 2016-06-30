program test
use jacobian
use global_parameters, only :PI
use binenergy
use binvel
use basic_functions
IMPLICIT NONE
REAL(KIND=8) :: x,y,z,h
external ERROR
REAL,DIMENSION(3) :: c_initial,c_final,wbin,fgh
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

call bin_velocity(fgh,x,c_initial,c_final,wbin)
call bin_energy(h,x,c_initial,c_final,wbin)
call bin_energy(y,x*1.00001,c_initial,c_final,wbin)
print *,wbin

print *,"The value of BINVEL is :",fgh
print *,"The value of BINENERGY is :",h
print *,"The value of the Jacobian is :",jac1D(x,c_initial,c_final,wbin)
print *,"The value of numerical Jac is:",(y-h)/(0.00001*x)
end program test
