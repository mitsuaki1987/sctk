!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_clock
  !
  IMPLICIT NONE
  !
CONTAINS
SUBROUTINE sctk_print_clock()
  !
  USE io_global, ONLY : stdout
  USE input_parameters, ONLY : calculation
  USE uspp, ONLY : okvan
  !
  USE sctk_val, ONLY : laddxc, lsf
  !
  IMPLICIT NONE
  !
  WRITE(stdout,'(/,5x,"#####  Clock summary  #####",/)')
  !
  IF(TRIM(calculation) == "kel") THEN
     !
     CALL print_clock("make_scrn")
     CALL print_clock("invert")
     CALL print_clock("make_kel")
     IF(laddxc>0) CALL print_clock("apply_xc")
     IF(lsf>0) CALL print_clock("lambda_sf")
     !
     WRITE(stdout,'(/,5x,"#####  In make_scrn #####",/)')
     CALL print_clock("fermi_factor")
     CALL print_clock("zgemm(make_scrn)")
     !
     WRITE(stdout,'(/,5x,"#####  In make_kel #####",/)')
     CALL print_clock("zgemm(make_kel)")
     !
     WRITE(stdout,'(/,5x,"#####  In make_scrn and make_kel #####",/)')
     CALL print_clock("wfc*wfc=rho")
     CALL print_clock("fft[rho]")
     !
     IF(okvan) THEN
        WRITE(stdout,'(/,5x,"#####  In wfc*wfc=rho #####",/)')
        CALL print_clock("addusxx")
     END IF
     !
  ELSE ! calculation = scdft, lambda_mu_k, qpdos, deltaf, ultrasonic
     !
     CALL print_clock("read_elph")
     CALL print_clock("read_Coulomb")
     !
     IF(TRIM(calculation) == "scdft") THEN
        !
        CALL print_clock("expand_g_v")
        CALL print_clock("make_Z")
        CALL print_clock("make_effint")
        CALL print_clock("broyden_gapeq")
        !
     ELSE ! calculation = lambda_mu_k, qpdos, deltaf, ultrasonic
        !
        CALL print_clock("expand_g_v_f")
        !
        IF(TRIM(calculation) == "lambda_mu_k") THEN
           !
           CALL print_clock("make_lambda_mu_f")
           !
        ELSE ! calculation = qpdos, deltaf, ultrasonic
           !
           IF(TRIM(calculation) == "qpdos") THEN
              !
              CALL print_clock("make_Z_qpdos")
              CALL print_clock("gapeq_rhs_qpdos")
              !
           ELSE ! deltaf, ultrasonic
              !
              CALL print_clock("make_Z_f")
              CALL print_clock("gapeq_rhs_f")
              !
           END IF
        END IF
     END IF
  END IF
  !
END SUBROUTINE sctk_print_clock
!
END MODULE sctk_clock
