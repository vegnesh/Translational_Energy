program onevar
  use global_parameters , only : K_B, ATOMMASS
  use binenergy
  use binvel
  use jacobian
  IMPLICIT NONE
  REAL :: beta, energy, resultv,jacval
  REAL, DIMENSION(3) :: c_f,c_i,w,velbin
  REAL :: T
  INTEGER,PARAMETER :: nbins = 2
  INTEGER :: NITER,I
  REAL :: IGUESS,NGUESS
  REAL,PARAMETER :: MASS = 4.002602*ATOMMASS
  w= (/0.D0, 0.D0, 0.D0 /) 
  c_i =  0*(/-1000.D0, -1000.D0, -1000.D0 /) 
  c_f =  (/10.D0, 10.D0, 10.D0/)
  ! beta at T  = 1000  = 1
  write(*, '(A)', ADVANCE = "NO") "Enter the value of beta:  "
  read(*,*) beta
  write(*, '(A)', ADVANCE = "NO") "Enter the value of IGUESS:  "
  read(*,*) IGUESS
  call bin_energy(energy,beta,c_i,c_f,w)
  call bin_velocity(velbin,beta,c_i,c_f,w)
  print *,velbin,energy,c_i,c_f
  !IGUESS = 0.2
  do i=1,30
    call bin_energy(resultv,IGUESS,c_i,c_f,w)
    jacval = jacenergy_beta(IGUESS,c_i,c_f,w)
    NGUESS = IGUESS - (resultv - energy)/jacval
    IGUESS = NGUESS;
    !print *, "True", beta, "ENERGY",ENERGY,"GUESS",IGUESS,"RESULT",RESULTv,"Jac", jacval
  end do 
print *, "True", beta, "ENERGY",ENERGY,"GUESS",IGUESS,"RESULT",RESULTv,"Jac", jacval


  end program onevar
!
!  write(*, '(A)', ADVANCE = "NO") "Enter the value of beta:  "
!  read(*,*) beta
!  write(*, '(A)', ADVANCE = "NO") "Enter the value of Min Velocity:  "
!  read(*,*) c_i
!  write(*, '(A)', ADVANCE = "NO") "Enter the value of Max Velocity:  "
!  read(*,*) c_f
!  
!   
!    print *, "The value of the Jacobian  is " , jacenergy_beta(beta,c_i,c_f)
!    call bin_energy(resultv,beta,c_i,c_f)
!    call bin_energy(energy,beta*1.000001,c_i,c_f)
!    print *, "Numerical Jacobian" , (energy - resultv)/(0.000001*beta) 
!    print *, "The Boltzmann Constant is ", K_B

