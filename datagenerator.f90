program datagen
  use global_parameters , only : K_B, ATOMMASS
  use binenergy
  use jacobian
  IMPLICIT NONE
  REAL :: beta, energy, resultv,jacval
  REAL, DIMENSION(3) :: c_f,c_i,w,cilimit,cflimit,temp
  REAL :: width
  REAL,ALLOCATABLE,DIMENSION(:) :: ZerotoOne, OnetoN
  INTEGER :: NPOINTS,I,NBINS,statusval,as1,as2,J,K
  CHARACTER(len=40) :: Datafile, ichars,jchar, filename

  cilimit = (/-1.D14, -1.D14, -1.D14/)  
  cflimit = (/-3.D2, -3.D2, -3.D2/)  
  w = (/0.D0, 0.D0, 0.D0/)
  
!******************BETA VALUES GENERATION**************************    
  write(*, '(A)', ADVANCE = "NO") "Enter the number of values:  "
  read(*,*) NPOINTS
  ALLOCATE( zerotoone(1:(NPOINTS/10)),STAT = statusval )
  if(statusval .NE. 0) stop "Problem with zerotoone"
  ALLOCATE( OnetoN(1:(NPOINTS-NPOINTS/10)),STAT = statusval )
  if(statusval .NE. 0) stop "Problem with OnetoN"
  
  as1 = size(zerotoone)
  as2 = size(OnetoN,1)

  width = 1.0/as1
  do i=1,as1
    zerotoone(i) = (i)*width
  end do
  width = 15.0/as2 
  do i = 1,as2
    OnetoN(i) = exp((i)*width)
  end do

  
!********************DATA GENERATION******************************** 
  Datafile = "BINS_BIN"
  NBINS = 10 !!!! DATA GENERATED FOR 10 BINS

  do i=1,NBINS
    c_i = cilimit
    c_f = cflimit
    if (i==1) c_f = (-1)*cilimit
    if (i==2) c_f = 0*(-cilimit)

    do j = 1,i
    !  print *,"NBINS:",i,"BIN:",j,"C_i:",c_i,"c_f:",c_f   
!********************File handling********************************** 
      write(ichars,'(i2)') i
      write(jchar,'(i2)') j
      ichars = adjustl(ichars)
      jchar = adjustl(jchar)
      filename="DATA/"//trim(ichars)//trim(Datafile)//trim(jchar)//".dat"
      print *, filename
      open(2,file = filename)
      do k=1,as1
        call bin_energy(energy,zerotoone(k),c_i,c_f,w)
        write(2,'(f16.8,f16.8)') zerotoone(k), energy
      end do
      do k=1,as2
        call bin_energy(energy,OnetoN(k),c_i,c_f,w)
        write(2,'(f16.8,f16.8)') OnetoN(k), energy
      end do
      close(2)
!*******************************************************************    
      temp = c_i
      c_i = c_f
      if (mod(i,2)==0) then
         if (j<i/2 .and. j+1 .ne. i/2) then
           c_f = c_f/10
         elseif (j+1 == i/2) then
           c_f = 0*(-c_f)
         elseif (j==i/2) then
           c_f = -temp
         elseif (j==i-1) then
           c_f = -cilimit
         else
           c_f = c_f*10
         endif
       else
         if (j<i/2) then 
           c_f = c_f/10
         elseif (j==i/2) then
           c_f = -c_f
         elseif (j==i-1) then 
           c_f = -cilimit
         else 
           c_f = c_f*10
         endif  
       endif

    end do
  end do

!******************************************************************* 
 
  DEALLOCATE(zerotoone)
  DEALLOCATE(OnetoN)  
end program datagen 
