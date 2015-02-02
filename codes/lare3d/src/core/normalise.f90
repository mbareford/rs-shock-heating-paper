!**************************************************************
! This module contains the conversion factors between normalised
! and internal code units
!**************************************************************

MODULE normalise
  USE constants
  USE shared_data
  IMPLICIT NONE
  
CONTAINS

  
  SUBROUTINE normalise_transport
  
    ! Normalise tbar, rbar and etabar for including Cowling resistivity and neutrals
    tbar = tbar/Tp0
    rbar = rbar*Dn0/(Tp0**1.5_num)  
    etabar = etabar/eta0
 
    ! Normalise ionise_pot 
    ionise_pot = ionise_pot_si/(En0*m0)
    IF (eos_number /= EOS_ION) ionise_pot = 0.0_num
         
    ! normalise tr required for get_neutral etc.
    tr = tr/Tp0            
  
  END SUBROUTINE normalise_transport
  


END MODULE normalise
