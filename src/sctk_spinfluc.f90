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
  USE mp_world, ONLY : mpime, world_comm
  !
  USE sctk_dmuxc, ONLY : apply_xc_spin
  USE sctk_invert, ONLY : invert, hermite
  USE noncollin_module, ONLY : npol
  USE mp, ONLY : mp_bcast
  !
  USE sctk_val, ONLY : nmf, wscr, ngv, ngv0, ngv1
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ipol
  REAL(dp) :: stoner(2:2*npol)
  LOGICAL :: approx
  !
  approx = .false.
  !
  CALL start_clock("lambda_sf")
  !
  CALL apply_xc_spin()
  !
  ! Stoner factor
  !
  IF(mpime == 0) THEN
     stoner(2:2*npol) = REAL(wscr(1,1,0,2:2*npol), dp) 
     WRITE(*,'(/,9x,"Stoner factor:")')
     WRITE(*,'(/,11x,3(e12.5,2x))') stoner(2:2*npol)
  END IF
  CALL mp_bcast(stoner, 0, world_comm)
  !
  CALL hermite()
  !
  IF(approx) THEN
     !
     DO ipol = 2, 2*npol
        wscr(1:ngv,ngv0:ngv1,0:nmf,ipol) = wscr(1:ngv,ngv0:ngv1,0:nmf,ipol) &
        &                                / (1.0_dp - stoner(ipol))
     END DO
     !
  ELSE
     !
     wscr(1:ngv,ngv0:ngv1,0:nmf,2:2*npol) = - wscr(1:ngv,ngv0:ngv1,0:nmf,2:2*npol)
     DO ig = ngv0, ngv1
        wscr(ig,ig,0:nmf,2:2*npol) = 1.0_dp + wscr(ig,ig,0:nmf,2:2*npol)
     END DO
     !
     CALL apply_xc_spin()
     !
     DO ipol = 2, 2*npol
        CALL invert(wscr(1:ngv,ngv0:ngv1,0:nmf,ipol))
     END DO
     !
     CALL apply_xc_spin()
     !
     DO ig = ngv0, ngv1
        wscr(ig,ig,0:nmf,2:2*npol) = wscr(ig,ig,0:nmf,2:2*npol) - 1.0_dp
     END DO
     !
     CALL hermite()
     !
  END IF
  !
  CALL apply_xc_spin()
  !
  CALL hermite()
  !
  IF(npol == 1) THEN
     wscr(1:ngv,ngv0:ngv1,0:nmf,2) = - 3.0_dp*wscr(1:ngv,ngv0:ngv1,0:nmf,2)
  ELSE
     wscr(1:ngv,ngv0:ngv1,0:nmf,2:2*npol) = - wscr(1:ngv,ngv0:ngv1,0:nmf,2:2*npol)
  END IF
  !
  CALL stop_clock("lambda_sf")
  !
END SUBROUTINE lambda_sf
!
END MODULE sctk_spinfluc
