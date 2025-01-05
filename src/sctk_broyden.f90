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
!>
!>
SUBROUTINE output_Zcp(iter,Zcp)
  !
  USE kinds, ONLY : dp
  USE sctk_val, ONLY: nmf, mf
  !
  INTEGER,INTENT(IN) :: iter
  REAL(DP),INTENT(IN) :: Zcp(nmf,3)
  !
  INTEGER :: imf
  !
  DO imf = 1, nmf
    WRITE(200,*) iter, mf(imf), Zcp(imf,1), Zcp(imf,2), Zcp(imf,3)
  END DO
  WRITE(200,*) ""
  !
END SUBROUTINE output_Zcp
!>
!>
!>
SUBROUTINE compute_dabs(lout,delta,dabs)
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE constants, ONLY : RYTOEV
  USE sctk_val, ONLY : dk, emin, ngap, ngap1, ngap2, xi
  !
  LOGICAL,INTENT(IN) :: lout
  REAL(DP),INTENT(INOUT) :: delta(ngap,2)
  REAL(dp),INTENT(OUT) :: dabs
  !
  LOGICAL :: lfermi(ngap)
  REAL(dp) :: ave(2), dave(2)
  !
  ave(1:2) = 0.0_dp
  dave(1:2) = 0.0_dp
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
  IF(lout) THEN
    WRITE(stdout,'(9x," Delta [meV] : ",2(e12.5,2x))') dave(1:2)
  ELSE
    delta(1:ngap,1) = delta(1:ngap,1) * SIGN(1.0_dp, dave(1))
    delta(1:ngap,2) = delta(1:ngap,2) * SIGN(1.0_dp, dave(2))
  END IF
  !
END SUBROUTINE
!>
!> Solve gap equation with the Broyden method
!>
SUBROUTINE broyden_gapeq(ndim,vec)
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE control_flags, ONLY : niter, tr2
  USE input_parameters, ONLY : calculation
  !
  USE sctk_gapeq_rhs, ONLY : gapeq_rhs
  USE sctk_eliashberg_sub, ONLY : eliashberg_rhs
  !
  REAL(dp),INTENT(OUT) :: vec(ndim)
  INTEGER,INTENT(IN) :: ndim
  !
  INTEGER :: itr, jtr
  REAL(dp) :: out_in_norm, mix = 0.2_dp, d_vec(ndim), out_in(ndim), out_in_old(ndim), d_out_in(ndim), &
  &          jacob1(ndim,2:niter), jacob2(ndim,2:niter)
  REAL(dp) :: dabs
  REAL(dp),EXTERNAL :: ddot
  !
  CALL start_clock("broyden_gapeq")
  !
  DO itr = 1, niter
    !
    WRITE(stdout,'(/,7x,"Iteration ",i0,/)') itr
    !
    IF(calculation == "scdft" .OR. calculation == "scdft_tc") THEN
      CALL gapeq_rhs(vec,out_in)
    ELSE
      CALL eliashberg_rhs(vec,out_in)
    END IF
    !
    out_in_norm = ddot(ndim, out_in, 1, out_in, 1)
    out_in_norm = SQRT(out_in_norm) / REAL(ndim, dp)
    !
    WRITE(stdout,'(9x,"Residual[Ry] : ",e12.5)') out_in_norm
    IF(calculation == "scdft" .OR. calculation == "scdft_tc") THEN
      CALL compute_dabs(.TRUE., vec, dabs)
    ELSE
      CALL output_Zcp(itr, vec)
    END IF
    !
    IF(out_in_norm < tr2) EXIT
    !
    IF(itr > 1) THEN
      !
      d_out_in(1:ndim) = out_in(1:ndim) - out_in_old(1:ndim)
      jacob1(1:ndim,itr) = mix * d_out_in(1:ndim) + d_vec(1:ndim)
      DO jtr = 2, itr - 1
        jacob1(1:ndim,itr) = jacob1(1:ndim,itr) - jacob1(1:ndim,jtr) &
        &          * ddot(ndim, jacob2(1:ndim,jtr), 1, d_out_in(1:ndim), 1)
      END DO
      jacob2(1:ndim,itr) = d_out_in(1:ndim) / ddot(ndim, d_out_in(1:ndim), 1, d_out_in(1:ndim), 1)
      !
    END IF
    !
    d_vec(1:ndim) = mix * out_in(1:ndim)
    DO jtr = 2, itr
       d_vec(1:ndim) = d_vec(1:ndim) - jacob1(1:ndim,jtr) &
       &        * ddot(ndim, jacob2(1:ndim,jtr), 1, out_in(1:ndim), 1)
    END DO
    !
    out_in_old(1:ndim) = out_in(1:ndim)
    vec(1:ndim) = vec(1:ndim) + d_vec(1:ndim)
    !
  END DO ! itr
  !
  IF(itr >= niter) THEN
    WRITE(stdout,'(/,7x,"Not converged ! res = ",e12.5,/)') out_in_norm
  ELSE
    WRITE(stdout,'(/,7x,"Converged ! iter = ",i0,/)') itr
  END IF
  !
  CALL stop_clock("broyden_gapeq")
  !
END SUBROUTINE broyden_gapeq
!
END MODULE sctk_broyden
