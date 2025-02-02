!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_spinfluc
  !
  IMPLICIT NONE
  !
CONTAINS
!
SUBROUTINE lambda_sf()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime
  !
  USE sctk_dmuxc, ONLY : apply_xc_spin
  USE sctk_invert, ONLY : invert, hermite, eigmax_power
  USE noncollin_module, ONLY : npol
  USE mp, ONLY : mp_bcast
  !
  USE sctk_val, ONLY : nmf, wscr, ngv, ngv0, ngv1
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ipol, iter(2:2*npol)
  REAL(dp) :: stoner(2:2*npol)
  !
  CALL start_clock("lambda_sf")
  !
  CALL apply_xc_spin()
  !
  ! Stoner factor
  !
  IF(mpime == 0) THEN
     stoner(2:2*npol) = REAL(wscr(1,1,0,2:2*npol), dp) 
     WRITE(*,'(/,9x,"Stoner criteria:")')
     WRITE(*,'(11x,"(I_{xc} Pi_0)_{00} : ",3(e12.5,2x))') stoner(2:2*npol)
  END IF
  DO ipol = 2, 2*npol
     CALL eigmax_power(wscr(1:ngv,ngv0:ngv1,0,ipol), stoner(ipol), iter(ipol))
  END DO
  IF(mpime == 0) THEN
     WRITE(*,'(11x,"EigMax(I_{xc} Pi_0) : ",3(e12.5,2x),/)') stoner(2:2*npol)
     WRITE(*,'(13x,"Iteration for power method : ",3(i6),/)') iter(2:2*npol)
  END IF
  !
  DO ig = ngv0, ngv1
     wscr(ig,ig,0:nmf,2:2*npol) = wscr(ig,ig,0:nmf,2:2*npol) - 1.0_dp
  END DO
  !
  DO ipol = 2, 2*npol
     CALL invert(wscr(1:ngv,ngv0:ngv1,0:nmf,ipol))
  END DO
  !
  DO ig = ngv0, ngv1
     wscr(ig,ig,0:nmf,2:2*npol) = wscr(ig,ig,0:nmf,2:2*npol) + 1.0_dp
  END DO
  !
  CALL hermite()
  !
  CALL apply_xc_spin()
  !
  IF(npol == 1) wscr(1:ngv,ngv0:ngv1,0:nmf,2) = 3.0_dp*wscr(1:ngv,ngv0:ngv1,0:nmf,2)
  !
  CALL stop_clock("lambda_sf")
  !
END SUBROUTINE lambda_sf
!
END MODULE sctk_spinfluc
