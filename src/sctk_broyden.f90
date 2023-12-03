!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_broyden
  !
  IMPLICIT NONE
  !
CONTAINS
!>
!> Solve gap equation with the Broyden method
!>
SUBROUTINE broyden_gapeq(lout_delta,dabs)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime
  USE io_global, ONLY : stdout
  USE control_flags, ONLY : niter, tr2
  USE sctk_val, ONLY : delta, dk, emin, ngap, ngap1, ngap2, xi
  USE constants, ONLY : RYTOEV
  USE io_files, ONLY : prefix, tmp_dir
  !
  USE sctk_gapeq_rhs, ONLY : gapeq_rhs
  USE sctk_io_delta, ONLY : out_delta
  !
  LOGICAL,INTENT(IN) :: lout_delta
  REAL(dp),INTENT(OUT) :: dabs
  !
  INTEGER :: itr, jtr
  REAL(dp) :: res, alpha = 0.2_dp, dd(ngap,2), rhs(ngap,2), rhs0(ngap,2), drhs(ngap,2), &
  &          jacob1(ngap,2,niter), jacob2(ngap,2,niter), ave(2), dave(2)
  LOGICAL :: lfermi(ngap)
  REAL(dp),EXTERNAL :: ddot
  !
  CALL start_clock("broyden_gapeq")
  !
  itr = 0
  WRITE(stdout,'(/,7x,"Iteration ",i0,/)') itr
  !
  CALL gapeq_rhs(rhs)
  res = ddot(ngap * 2, rhs, 1, rhs, 1)
  res = SQRT(res) / REAL(ngap * 2, dp)
  !
  ! $$$$$  Information of conversience  $$$$$
  !
  ave(1:2) = 0.0_dp
  dave(1:2) = 0.0_dp
  !
  lfermi(1:ngap1) = ABS(xi(1:ngap1,1)) < emin
  dave(1) = SUM(delta(1:ngap1,1) * dk(1:ngap1,1), lfermi(1:ngap1))
  ave(1)  = SUM(                 dk(1:ngap1,1), lfermi(1:ngap1))
  !
  lfermi(1:ngap2) = ABS(xi(1:ngap2,2)) < emin
  dave(2) = SUM(delta(1:ngap2,2) * dk(1:ngap2,2), lfermi(1:ngap2))
  ave(2)  = SUM(                 dk(1:ngap2,2), lfermi(1:ngap2))
  !
  dave(1:2) = dave(1:2) / ave(1:2) * RYTOEV * 1.0e3_dp
  !
  WRITE(stdout,'(9x,"Residual[Ry] : ",e12.5)') res
  WRITE(stdout,'(9x," Delta [meV] : ",2(e12.5,2x))') dave(1:2)
  !
  ! $$$$$  End information of conversience  $$$$$
  !
  IF(res < tr2) GOTO 5
  !
  dd(1:ngap,1:2) = - alpha * rhs(1:ngap,1:2)
  !
  DO itr = 1, niter
     !
     WRITE(stdout,'(/,7x,"Iteration ",i0,/)') itr
     !
     delta(1:ngap,1:2) = delta(1:ngap,1:2) + dd(1:ngap,1:2)
     !
     rhs0(1:ngap,1:2) = rhs(1:ngap,1:2)
     CALL gapeq_rhs(rhs)
     res = ddot(ngap * 2, rhs, 1, rhs, 1)
     res = SQRT(res) / REAL(ngap, dp)
     !
     ! $$$$$  Information of conversience  $$$$$
     !
     ave(1:2) = 0.0_dp
     dave(1:2) = 0.0_dp
     !
     lfermi(1:ngap1) = ABS(xi(1:ngap1,1)) < emin
     dave(1) = SUM(delta(1:ngap1,1) * dk(1:ngap1,1), lfermi(1:ngap1))
     ave(1)  = SUM(                 dk(1:ngap1,1), lfermi(1:ngap1))
     !
     lfermi(1:ngap2) = ABS(xi(1:ngap2,2)) < emin
     dave(2) = SUM(delta(1:ngap2,2) * dk(1:ngap2,2), lfermi(1:ngap2))
     ave(2)  = SUM(                 dk(1:ngap2,2), lfermi(1:ngap2))
     !
     dave(1:2) = dave(1:2) / ave(1:2) * RYTOEV * 1.0e3_dp
     !
     WRITE(stdout,'(9x,"Residual[Ry] : ",e12.5)') res
     WRITE(stdout,'(9x," Delta [meV] : ",2(e12.5,2x))') dave(1:2)
     !
     ! $$$$$  End information of conversience  $$$$$
     !
     IF(res < tr2) THEN
        !       
        delta(1:ngap,1) = delta(1:ngap,1) * SIGN(1.0_dp, dave(1))
        delta(1:ngap,2) = delta(1:ngap,2) * SIGN(1.0_dp, dave(2))
        !
        GOTO 5
        !
     END IF
     !
     ! Update Jacobian with drhs
     !
     drhs(1:ngap,1:2) = rhs(1:ngap,1:2) - rhs0(1:ngap,1:2)
     !
     jacob1(1:ngap,1:2,itr) = - alpha * drhs(1:ngap,1:2)
     DO jtr = 1, itr - 1
        jacob1(1:ngap,1:2,itr) = jacob1(1:ngap,1:2,itr) - jacob1(1:ngap,1:2,jtr) &
        &          * ddot(ngap * 2, jacob2(1:ngap,1:2,jtr), 1, drhs(1:ngap,1:2), 1)
     END DO
     jacob1(1:ngap,1:2,itr) = dd(1:ngap,1:2) + jacob1(1:ngap,1:2,itr)
     jacob2(1:ngap,1:2,itr) = drhs(1:ngap,1:2) / ddot(ngap*2, drhs(1:ngap,1:2), 1, drhs(1:ngap,1:2), 1)
     !
     ! Compute dd with new Jacobian & rhs
     !
     dd(1:ngap,1:2) = - alpha * rhs(1:ngap,1:2)
     DO jtr = 1, itr
        dd(1:ngap,1:2) = dd(1:ngap,1:2) - jacob1(1:ngap,1:2,jtr) &
        &        * ddot(ngap*2, jacob2(1:ngap,1:2,jtr), 1, rhs(1:ngap,1:2), 1)
     END DO
     !
  END DO ! itr
  !
  lfermi(1:ngap1) = ABS(xi(1:ngap1,1)) < emin
  dave(1) = SUM(ABS(delta(1:ngap1,1)) * dk(1:ngap1,1), lfermi(1:ngap1))
  ave(1)  = SUM(                        dk(1:ngap1,1), lfermi(1:ngap1))
  !
  lfermi(1:ngap2) = ABS(xi(1:ngap2,2)) < emin
  dave(2) = SUM(ABS(delta(1:ngap2,2)) * dk(1:ngap2,2), lfermi(1:ngap2))
  ave(2)  = SUM(                        dk(1:ngap2,2), lfermi(1:ngap2))
  !
  dave(1:2) = dave(1:2) / ave(1:2) * RYTOEV * 1.0e3_dp
  dabs = 0.5_dp * SUM(dave(1:2))
  !
  IF(mpime == 0) THEN
     IF(lout_delta) CALL out_delta(TRIM(tmp_dir) // TRIM(prefix) // ".scgap")
  END IF
  !
  WRITE(stdout,'(/,7x,"Not converged ! res = ",e12.5,/)') res
  RETURN
  !
5 CONTINUE
  !
  lfermi(1:ngap1) = ABS(xi(1:ngap1,1)) < emin
  dave(1) = SUM(ABS(delta(1:ngap1,1)) * dk(1:ngap1,1), lfermi(1:ngap1))
  ave(1)  = SUM(                        dk(1:ngap1,1), lfermi(1:ngap1))
  !
  lfermi(1:ngap2) = ABS(xi(1:ngap2,2)) < emin
  dave(2) = SUM(ABS(delta(1:ngap2,2)) * dk(1:ngap2,2), lfermi(1:ngap2))
  ave(2)  = SUM(                        dk(1:ngap2,2), lfermi(1:ngap2))
  !
  dave(1:2) = dave(1:2) / ave(1:2) * RYTOEV * 1.0e3_dp
  dabs = 0.5_dp * SUM(dave(1:2))
  !
  IF(mpime == 0) THEN
     IF(lout_delta) CALL out_delta(TRIM(tmp_dir) // TRIM(prefix) // ".scgap")
  END IF
  !
  WRITE(stdout,'(/,7x,"Converged ! iter = ",i0,/)') itr
  !
  CALL stop_clock("broyden_gapeq")
  !
END SUBROUTINE broyden_gapeq
!
END MODULE sctk_broyden
