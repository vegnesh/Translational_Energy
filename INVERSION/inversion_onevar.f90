program onevar
  use global_parameters , only : K_B, ATOMMASS
  use binenergy
  use binvel
  use jacobian

  IMPLICIT NONE
  REAL :: beta, energy, resultv,jacval
  REAL, DIMENSION(1:3) :: c_f,c_i,w,velbin
  REAL :: T
  INTEGER,PARAMETER :: nbins = 2
  INTEGER :: NITER,I
  REAL :: IGUESS,NGUESS
  REAL,PARAMETER :: MASS = 4.002602*ATOMMASS
  CHARACTER :: choice = 'y'

  do while (choice == 'y')
  w= (/0.D0, 0.D0, 0.D0 /) 
  c_i =  10.D0*(/-1000.D0, -1000.D0, -1000.D0 /) 
  c_f =  10.D3*(/10.D0, 10.D0, 10.D0/)
  ! beta at T  = 1000  = 1
  write(*, '(A)', ADVANCE = "NO") "Enter the value of beta:  "
  read(*,*) beta
  write(*, '(A)', ADVANCE = "NO") "Enter the value of Min Velocity:  "
  read(*,*) c_i(1)
  write(*, '(A)', ADVANCE = "NO") "Enter the value of Max Velocity:  "
  read(*,*) c_f(1)
  write(*, '(A)', ADVANCE = "NO") "Enter the value of IGUESS:  "
  read(*,*) IGUESS
  

  call bin_energy(energy,beta,c_i,c_f,w)
  call bin_velocity(velbin,beta,c_i,c_f,w)
  
  print *,energy
  !IGUESS = 0.2
  do i=1,40
    call bin_energy(resultv,IGUESS,c_i,c_f,w)
    jacval = jacenergy_beta(IGUESS,c_i,c_f,w)
    NGUESS = IGUESS - (resultv - energy)/jacval
    IGUESS = NGUESS;
    !print *, "True", beta, "ENERGY",ENERGY,"GUESS",IGUESS,"RESULT",RESULTv,"Jac", jacval
  end do 
  print *, "True", beta, "ENERGY",ENERGY,"GUESS",IGUESS,"RESULT",RESULTv,"Jac", jacval
  write(*, '(A)', ADVANCE = "NO") "Continue ?y?n? :  "
  read(*,*) choice
  
  end do

  end program onevar

