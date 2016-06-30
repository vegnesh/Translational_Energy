program onevar
  use global_parameters , only : K_B
  use binenergy
  IMPLICIT NONE
  REAL,EXTERNAL :: jac1D
  REAL :: beta, energy, resultv
  REAL, DIMENSION(3) :: c_f,c_i
  REAL :: T
  INTEGER,PARAMETER :: nbins = 2
  INTEGER :: NITER,I
  REAL :: IGUESS,NGUESS
  
  write(*, '(A)', ADVANCE = "NO") "Enter the value of beta:  "
  read(*,*) beta
  write(*, '(A)', ADVANCE = "NO") "Enter the value of Min Velocity:  "
  read(*,*) c_i
  write(*, '(A)', ADVANCE = "NO") "Enter the value of Max Velocity:  "
  read(*,*) c_f
  
  ! ONE BIN SIMULATION
  !c_i(1:3) = (/-1.D15, -1.D15, -1.D15/)
  !c_f(1:3) = (/1.D15, 1.D15, 1.D15/)
 
 print *, "The value of the Jacobian  is " , jac1D(beta,c_i,c_f)
    call bin_energy(resultv,beta,c_i,c_f)
    call bin_energy(energy,beta*1.000001,c_i,c_f)
    print *, "Numerical Jacobian" , (energy - resultv)/(0.000001*beta) 
    print *, "The Boltzmann Constant is ", K_B


 ! do I = 1,1
 !   T = 200.D0 + 50.D0 * I
 !   beta = 1/(K_B*T)
 !   IGUESS = 210.D0
 !   NITER = 0
 !   do while (NITER .lt. 100 )!.or. ABS(T*1.0-IGUESS)/IGUESS>1.E-6)
 !     call bin_energy(energy,(1.0/(K_B*IGUESS)),c_i,c_f)
 !     NGUESS = IGUESS - energy/jac1D(IGUESS,c_i,c_f)
 !     IGUESS = NGUESS
 !     NITER = NITER +1
 !     print *, IGUESS
 !    !   end do
 !   print *,"T:",T,"NITER:",NITER,"Sol:", IGUESS,"Act:", beta
 ! end do
end program onevar
