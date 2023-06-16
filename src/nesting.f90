!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE nesting_tetra
  !
  IMPLICIT NONE
  !
CONTAINS
!
SUBROUTINE nesting_delta1(nbnd_fs,nspin_lsda,nirr_k,irr_k,eig1,eig2,nesting_factor)
  !---------------------------------------------------------------------
  !
  ! This routine computed the weight for the double-delta function.
  !
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE klist,  ONLY: nkstot
  USE lsda_mod,   ONLY : nspin
  USE ktetra, ONLY : tetra, ntetra, nntetra, wlsm
  !
  INTEGER,INTENT(IN) :: nbnd_fs, nspin_lsda, nirr_k
  INTEGER,INTENT(IN) :: irr_k(3,nirr_k)
  REAL(dp),INTENT(IN) :: eig1(nbnd_fs, nkstot, nspin_lsda), eig2(nbnd_fs*nirr_k,nkstot, nspin_lsda)
  REAL(dp),INTENT(OUT) :: nesting_factor(nirr_k,nspin_lsda)
  !
  INTEGER :: ns, nt, ii, ibnd, itetra(4), start_t, last_t
  REAL(dp) :: e(4), a(4,4), V, tsmall(3,4)
  REAL(dp),ALLOCATABLE :: ei0(:,:), ej0(:,:), ej1(:,:), nesting_factor1(:)
  !
  nesting_factor(1:nirr_k,1:nspin_lsda) = 0.0_dp
  !
  CALL divide(world_comm, 6 * nkstot, start_t, last_t)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & PRIVATE(ns,nt,ibnd,ii,e,ei0,ej0,a,V,tsmall,ej1,itetra,nesting_factor1) &
  !$OMP & SHARED(nspin_lsda,nbnd_fs,nirr_k,start_t,last_t,eig1,eig2,tetra,wlsm,nntetra,irr_k) &
  !$OMP & REDUCTION(+:nesting_factor)
  !
  ALLOCATE(ei0(4,nbnd_fs), ej0(4,nbnd_fs*nirr_k), ej1(3,nbnd_fs*nirr_k), nesting_factor1(nirr_k))
  !
  DO ns = 1, nspin_lsda
     !
     !$OMP DO
     DO nt = start_t, last_t
        !
        ei0(1:4, 1:nbnd_fs) = 0.0_dp
        ej0(1:4, 1:nbnd_fs*nirr_k) = 0.0_dp
        DO ii = 1, nntetra
           !
           DO ibnd = 1, nbnd_fs
              ei0(1:4, ibnd) = ei0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * eig1(ibnd, tetra(ii, nt), ns)
           END DO
           !
           DO ibnd = 1, nbnd_fs*nirr_k
              ej0(1:4, ibnd) = ej0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * eig2(ibnd, tetra(ii, nt), ns)
           END DO
           !
        END DO
        !
        DO ibnd = 1, nbnd_fs
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(1:4,ii) = (0.0_dp - e(ii)) / (e(1:4) - e(ii))
           END DO
           !
           IF(e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) THEN
              !
              ! A
              !
              !V = 3.0_dp * a(2,1) * a(3,1) * a(4,1) / (0.0_dp - e(1))
              V = 3.0_dp * a(2,1) * a(3,1)           / (e(4) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              !
              ej1(1:3,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
              !
              CALL nesting_delta2(nbnd_fs,nirr_k,irr_k,ej1,nesting_factor1)
              !
              nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
              !
           ELSE IF( e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) THEN
              !
              ! B - 1
              !
              !V = 3.0_dp * a(3,1) * a(4,1) * a(2,4) / (0.0_dp - e(1))
              V = 3.0_dp           * a(4,1) * a(2,4) / (e(3) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              !
              ej1(1:3,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
              !
              CALL nesting_delta2(nbnd_fs,nirr_k,irr_k,ej1,nesting_factor1)
              !
              nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
              !
              ! B - 2
              !
              !V = 3.0_dp * a(2,3) * a(3,1) * a(4,2) / (0.0_dp - e(1))
              V = 3.0_dp * a(2,3)           * a(4,2) / (e(3) - e(1))
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              !
              ej1(1:3,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
              !
              CALL nesting_delta2(nbnd_fs,nirr_k,irr_k,ej1,nesting_factor1)
              !
              nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
              !
           ELSE IF(e(3) < 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C
              !
              !V = 3.0_dp * a(1,4) * a(2,4) * a(3,4) / (e(4) - 0.0_dp)
              V = 3.0_dp * a(1,4) * a(2,4)           / (e(4) - e(3))
              !
              tsmall(1, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
              !
              ej1(1:3,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:3,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
              !
              CALL nesting_delta2(nbnd_fs,nirr_k,irr_k,ej1,nesting_factor1)
              !
              nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
              !
           END IF
           !
        END DO
        !
     END DO ! nt
     !$OMP END DO
     !
  END DO ! ns
  !
  DEALLOCATE(ei0, ej0, ej1, nesting_factor1)
  !
  !$OMP END PARALLEL
  !
  nesting_factor(1:nirr_k,1:nspin_lsda) = nesting_factor(1:nirr_k,1:nspin_lsda) / REAL(ntetra, dp)
  IF(nspin == 1) nesting_factor(1:nirr_k,1) = 2.0_dp * nesting_factor(1:nirr_k,1)
  !
  CALL mp_sum(nesting_factor, world_comm)
  !
END SUBROUTINE nesting_delta1
!
!-------------------------------------------------------------------------------
SUBROUTINE nesting_delta2(nbnd_fs,nirr_k,irr_k,ej0,nesting_factor)
  !-----------------------------------------------------------------------------
  !
  ! 2nd step of tetrahedra method.
  !
  USE kinds, ONLY : dp
  !
  INTEGER,INTENT(IN) :: nbnd_fs, nirr_k
  INTEGER,INTENT(IN) :: irr_k(3,nirr_k)
  REAL(dp),INTENT(IN) :: ej0(3,nbnd_fs,nirr_k)
  REAL(dp),INTENT(OUT) :: nesting_factor(nirr_k)
  !
  INTEGER :: ibnd, itetra(3), ii, ik
  REAL(dp) :: e(3), a(3,3), V
  !
  nesting_factor(1:nirr_k) = 0.0_dp
  !
  DO ik = 1, nirr_k
     !
     IF(ALL(irr_k(1:3,ik) == 0)) CYCLE
     !
     DO ibnd = 1, nbnd_fs
        !
        IF(MAXVAL(ABS(ej0(1:3,ibnd,ik))) < 1e-10_dp) &
        & CALL errore("nesting_delta2", "Nesting occurs.", ibnd)
        !
        itetra(1) = 0
        e(1:3) = ej0(1:3,ibnd,ik)
        call hpsort (3, e, itetra)
        !
        DO ii = 1, 3
           a(1:3,ii) = (0.0_dp - e(ii)) / (e(1:3) - e(ii))
        END DO
        !
        IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
           !
           !V = 2.0_dp * a(2,1) * a(3,1) / (0.0_dp - e(1))
           V = 2.0_dp * a(2,1)           / (e(3) - e(1))
           !
           nesting_factor(ik) = nesting_factor(ik) + V
           !
        ELSE IF((e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) .OR. (e(2) < 0.0_dp .AND. 0.0_dp <= e(3))) THEN
           !
           !V = 2.0_dp * a(1,3) * a(2,3) / (e(3) - 0.0_dp)
           V = 2.0_dp * a(1,3)           / (e(3) - e(2))
           !
           nesting_factor(ik) = nesting_factor(ik) + V
           !
        END IF
        !
     END DO ! ib
     !
  END DO
  !
END SUBROUTINE nesting_delta2
!
SUBROUTINE nesting_theta1(nbnd_fs,nspin_lsda,nirr_k,eig1,eig2,nesting_factor)
  !---------------------------------------------------------------------
  !
  ! This routine computed the weight for the double-delta function.
  !
  USE kinds, ONLY : dp
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE klist,  ONLY: nkstot
  USE lsda_mod,   ONLY : nspin
  USE ktetra, ONLY : tetra, ntetra, nntetra, wlsm
  !
  INTEGER,INTENT(IN) :: nbnd_fs, nspin_lsda, nirr_k
  REAL(dp),INTENT(IN) :: eig1(nbnd_fs, nkstot, nspin_lsda), eig2(nbnd_fs*nirr_k,nkstot, nspin_lsda)
  REAL(dp),INTENT(OUT) :: nesting_factor(nirr_k,nspin_lsda)
  !
  INTEGER :: ns, nt, ii, ibnd, itetra(4), start_t, last_t
  REAL(dp) :: e(4), a(4,4), thr = 1.0e-8_dp, V, tsmall(4,4), ei1(4)
  REAL(dp),ALLOCATABLE :: ei0(:,:), ej0(:,:), ej1(:,:), nesting_factor1(:)
  !
  nesting_factor(1:nirr_k,1:nspin_lsda) = 0.0_dp
  !
  CALL divide(world_comm, 6 * nkstot, start_t, last_t)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & PRIVATE(ns,nt,ibnd,ii,e,ei0,ej0,a,V,tsmall,ei1,ej1,itetra,nesting_factor1) &
  !$OMP & SHARED(nspin_lsda,nbnd_fs,nirr_k,start_t,last_t,eig1,eig2,tetra,wlsm,nntetra,thr) &
  !$OMP & REDUCTION(+:nesting_factor)
  !
  ALLOCATE(ei0(4,nbnd_fs), ej0(4,nbnd_fs*nirr_k), ej1(4,nbnd_fs*nirr_k), nesting_factor1(nirr_k))
  !
  DO ns = 1, nspin_lsda
     !
     !$OMP DO
     DO nt = start_t, last_t
        !
        ei0(1:4, 1:nbnd_fs) = 0.0_dp
        ej0(1:4, 1:nbnd_fs*nirr_k) = 0.0_dp
        DO ii = 1, nntetra
           !
           DO ibnd = 1, nbnd_fs
              ei0(1:4, ibnd) = ei0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * eig1(ibnd, tetra(ii, nt), ns)
           END DO
           !
           DO ibnd = 1, nbnd_fs*nirr_k
              ej0(1:4, ibnd) = ej0(1:4, ibnd) &
              &             + wlsm(1:4,ii) * eig2(ibnd, tetra(ii, nt), ns)
           END DO
           !
        END DO
        !
        DO ibnd = 1, nbnd_fs
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(1:4,ii) = (0.0_dp - e(ii)) / (e(1:4) - e(ii))
           END DO
           !
           IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(4, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
           ELSE IF( e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
           ELSE IF( e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd_fs*nirr_k) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd_fs*nirr_k))
                 !
                 CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
                 !
                 nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + V * nesting_factor1(1:nirr_k)
                 !
              END IF
              !
           ELSE IF( e(4) <= 0.0_dp ) THEN
              !
              ! D - 1
              !
              ei1(1:4) = e(1:4)
              ej1(1:4,1:nbnd_fs*nirr_k) = ej0(itetra(1:4), 1:nbnd_fs*nirr_k)
              !
              CALL nesting_theta2(nbnd_fs,nirr_k,ei1,ej1,nesting_factor1)
              !
              nesting_factor(1:nirr_k,ns) = nesting_factor(1:nirr_k,ns) + nesting_factor1(1:nirr_k)
              !
           END IF
           !
        END DO
        !
     END DO ! nt
     !$OMP END DO
     !
  END DO ! ns
  !
  DEALLOCATE(ei0, ej0, ej1, nesting_factor1)
  !
  !$OMP END PARALLEL
  !
  nesting_factor(1:nirr_k,1:nspin_lsda) = nesting_factor(1:nirr_k,1:nspin_lsda) / REAL(ntetra, dp)
  IF(nspin == 1) nesting_factor(1:nirr_k,1) = 2.0_dp * nesting_factor(1:nirr_k,1)
  !
  CALL mp_sum(nesting_factor, world_comm)
  !
END SUBROUTINE nesting_theta1
!
!-------------------------------------------------------------------------------
SUBROUTINE nesting_theta2(nbnd_fs,nirr_k,ei0,ej0,nesting_factor)
  !-----------------------------------------------------------------------------
  !
  ! 2nd step of tetrahedra method.
  !
  USE kinds, ONLY : dp
  !
  INTEGER,INTENT(IN) :: nbnd_fs, nirr_k
  REAL(dp),INTENT(IN) :: ei0(4), ej0(4,nbnd_fs,nirr_k)
  REAL(dp),INTENT(OUT) :: nesting_factor(nirr_k)
  !
  INTEGER :: ibnd, itetra(4), ii, ik
  REAL(dp) :: e(4), a(4,4), V, ei1(4), ej1(4), thr = 1.0e-8_dp, &
  &           tsmall(4,4), nesting_factor1
  !
  nesting_factor(1:nirr_k) = 0.0_dp
  !
  DO ik = 1, nirr_k
     DO ibnd = 1, nbnd_fs
        !
        itetra(1) = 0
        e(1:4) = ej0(1:4,ibnd,ik)
        call hpsort (4, e, itetra)
        !
        DO ii = 1, 4
           a(1:4,ii) = (0.0_dp - e(ii)) / (e(1:4) - e(ii))
        END DO
        !
        IF(0_dp <= e(1)) THEN
             !
             ! A - 1
             !
             ej1(1:4) = e(1:4)
             ei1(1:4) = ei0(itetra(1:4))
             !
             CALL nesting_lindhard(ei1,ej1,nesting_factor1)
             !
             nesting_factor(ik) = nesting_factor(ik) + nesting_factor1
             !
        ELSE IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
           !
           ! B - 1
           !
           V = a(1,2)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
              tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !       
           END IF
           !
           ! B - 2
           !
           V = a(1,3) * a(2,1)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !
           END IF
           !
           ! B - 3
           !
           V = a(1,4) * a(2,1) * a(3,1)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !       
           END IF
           !          
        ELSE IF((e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) .OR. (e(2) <= 0.0_dp .AND. 0.0_dp < e(3))) THEN
           !          
           ! C - 1
           !
           V = a(2,4) * a(1,4) * a(3,1)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !      
           END IF
           !
           ! C - 2
           !
           V = a(1,3) * a(2,3)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !
           END IF
           !
           ! C - 3
           ! 
           V = a(1,3) * a(2,4) * a(3,2)
           !
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
              tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !
           END IF
           !          
        ELSE IF((e(3) < 0.0_dp .AND. 0.0_dp <= e(4)) .OR. (e(3) <= 0.0_dp .AND. 0.0_dp < e(4))) THEN
           !
           ! D - 1
           !
           V = a(3,4) * a(2,4) * a(1,4) 
           !          
           IF(V > thr) THEN
              !
              tsmall(1, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
              tsmall(2, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
              tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
              tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
              !
              ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
              ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
              !
              CALL nesting_lindhard(ei1,ej1,nesting_factor1)
              !
              nesting_factor(ik) = nesting_factor(ik) + V * nesting_factor1
              !        
           END IF
           !
        END IF        !
     END DO ! ib
     !
  END DO
  !
END SUBROUTINE nesting_theta2
  !
SUBROUTINE nesting_lindhard(ei,ej,nesting_factor)
  !---------------------------------------------------------------
  !
  ! This routine compute 1 / (e_{k+q} - e_{k})
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei(4), ej(4)
  REAL(dp),INTENT(OUT) :: nesting_factor
  !
  INTEGER :: ii, itetra(4)
  REAL(dp) :: e(4), le(4), thr, thr2, w
  !
  w = 0.0_dp
  !
  itetra(1) = 0
  e(1:4) = ej(1:4) - ei(1:4)
  call hpsort (4, e, itetra)
  !
  thr = MAXVAL(e(1:4)) * 1e-3_dp
  thr2 = 1e-13_dp
  !
  DO ii = 1, 4
     IF(e(ii) < thr2) THEN
        IF(ii == 3) THEN
           WRITE(*,*) e(1:4)
           CALL errore("nesting_lindhard", "Nesting occurs.", 1)
        END IF
        le(ii) = 0.0_dp
        e(ii) = 0.0_dp
     ELSE
        le(ii) = LOG(e(ii))
     END IF
  END DO
  !
  IF(e(4) - e(3) < thr ) THEN
     IF(e(4) - e(2) < thr ) THEN
        IF(e(4) - e(1) < thr ) THEN
           !
           ! e(4) = e(3) = e(2) = e(1)
           !
           w = (1.0_dp)/e(4)
           !
        ELSE
           !
           ! e(4) = e(3) = e(2)
           !
           CALL nesting_1114(e(4),e(1),le(4),le(1),w)
           !
           IF(w < 0.0_dp) THEN
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(100e15.5)') w
              CALL errore("nesting_lindhard", "Case: e4=e3=e2", 1)
           END IF
           !
        END IF
     ELSE IF(e(2) - e(1) < thr ) THEN
        !
        ! e(4) = e(3), e(2) = e(1)
        !
        CALL nesting_1144(e(2),e(4),le(2),le(4),w)
        !
        IF(w < 0.0_dp) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w
           CALL errore("nesting_lindhard", "Case: e4=e3/=e2=e1", 1)
        END IF
        !
     ELSE
        !
        ! e(4) = e(3)
        !
        CALL nesting_1134(e(4),e(1),e(2),le(4),le(1),le(2),w)
        !
        IF(w < 0.0_dp) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w
           CALL errore("nesting_lindhard", "Case: e4=e3", 1)
        END IF
        !
     END IF
  ELSE IF(e(3) - e(2) < thr ) THEN
     IF(e(3) - e(1) < thr ) THEN
        !
        ! e(3) = e(2) = e(1)
        !
        CALL nesting_1114(e(3),e(4),le(3),le(4),w)
        !
        IF(w < 0.0_dp) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') ei(1:4)
           WRITE(*,'(100e15.5)') ej(1:4)
           WRITE(*,'(100e15.5)') w
           CALL errore("nesting_lindhard", "Case: e3=e2=e1", 1)
        END IF
        !
     ELSE
        !
        ! e(3) = e(2)
        !
        CALL nesting_1134(e(3),e(1),e(4),le(3),le(1),le(4),w)
        !
        IF(w < 0.0_dp) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w
           CALL errore("nesting_lindhard", "Case: e3=e2", 1)
        END IF
        !
     END IF
  ELSE IF(e(2) - e(1) < thr ) THEN
     !
     ! e(2) = e(1)
     !
     CALL nesting_1134(e(2),e(3),e(4),le(2),le(3),le(4),w)
     !
     IF(w < 0.0_dp) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w
        CALL errore("nesting_lindhard", "Case: e2=e1", 1)
     END IF
     !
  ELSE
     !
     ! Different each other.
     !
     CALL nesting_1234(e(1),e(2),e(3),e(4),le(1),le(2),le(3),le(4),w)
     !
     IF(w < 0.0_dp) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w
        CALL errore("nesting_lindhard", "Case: e1/=e2/=e3/=e4", 1)
     END IF
     !
  END IF
  !
  nesting_factor = w
  !
END SUBROUTINE nesting_lindhard
!
SUBROUTINE nesting_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4,w)
  !
  USE kinds, ONLY : dp
  !
  REAL(dp),INTENT(IN) :: g1, g2, g3, g4, lng1, lng2, lng3, lng4
  REAL(dp),INTENT(OUT) :: w
  !
  REAL(dp) :: w2, w3, w4
  !
  w2 = (g2*(lng2 - lng1))/(g2 - g1)
  w3 = (g3*(lng3 - lng1))/(g3 - g1)
  w4 = (g4*(lng4 - lng1))/(g4 - g1)
  w2 = (g2*(w2 - w3))/(g2 - g3)
  w4 = (g4*(w4 - w3))/(g4 - g3)
  w = (3.0_dp*(w4 - w2))/(g4 - g2)
  !
END SUBROUTINE nesting_1234
!
SUBROUTINE nesting_1134(g1,g3,g4,lng1,lng3,lng4,w)
  !
  USE kinds, ONLY : dp
  !
  REAL(dp),INTENT(IN) :: g1, g3, g4, lng1, lng3, lng4
  REAL(dp),INTENT(OUT) :: w
  !
  REAL(dp) :: w3, w4
  !
  w3 = (g3*(lng3 - lng1))/(g3 - g1)
  w4 = (g4*(lng4 - lng1))/(g4 - g1)
  w3 = (g3*(w3 - 1.0_dp))/(g3 - g1)
  w4 = (g4*(w4 - 1.0_dp))/(g4 - g1)
  w = (3.0_dp*(w4 - w3))/(g4 - g3)
  !
END SUBROUTINE nesting_1134
!
SUBROUTINE nesting_1114(g1,g4,lng1,lng4,w)
  !
  USE kinds, ONLY : dp
  !
  REAL(dp),INTENT(IN) :: g1, g4, lng1, lng4
  REAL(dp),INTENT(OUT) :: w
  !
  w = (g4*(lng4 - lng1))/(g4 - g1)
  w = (2.0_dp*g4*(w - 1.0_dp))/(g4 - g1)
  w = (3.0_dp*(w - 1.0_dp))/(2.0_dp*(g4 - g1))
  !
END SUBROUTINE nesting_1114
!
SUBROUTINE nesting_1144(g1,g4,lng1,lng4,w)
  !
  USE kinds, ONLY : dp
  !
  REAL(dp),INTENT(IN) :: g1, g4, lng1, lng4
  REAL(dp),INTENT(OUT) :: w
  !
  w = (g1*(lng4 - lng1))/(g4 - g1)
  w = (2.0_dp*g4*(1.0_dp - w))/(g4 - g1)
  w = (3.0_dp*(w - 1.0_dp))/(g4 - g1)
  !
END SUBROUTINE nesting_1144
!
END MODULE nesting_tetra
!
! Usage :
! $ nesting.x -in {pw.x input file}
! Then it generates nesting.frmsf (for nspin = 1, 4) or
! nesting.frmsf and nesting.frmsf (for nspin = 2)
!
!----------------------------------------------------------------------------
PROGRAM nesting
  !--------------------------------------------------------------------------
  !
  USE parameters,           ONLY : npk
  USE input_parameters,     ONLY : prefix, outdir
  USE io_files,             ONLY : prefix_ => prefix, tmp_dir
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE environment,          ONLY : environment_start, environment_end
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : et
  USE start_k,              ONLY : nk1, nk2, nk3
  USE cell_base,            ONLY : at, bg
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef, ef_up, ef_dw
  USE klist,                ONLY : nks, two_fermi_energies, xk, nkstot, wk
  USE fermisurfer_common,   ONLY : b_low, b_high, rotate_k_fs, write_fermisurfer
  USE ktetra,               ONLY : opt_tetra_init, tetra_type, opt_tetra_dos_t
  USE nesting_tetra,        ONLY : nesting_delta1, nesting_theta1
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i1, i2, i3, ns, nk, nbnd_fs, nirr_k, ikqv(3), s_dummy(3,3,48), t_rev_dummy(48)
  REAL(DP) :: ef1, ef2, dosf(2), edos
  INTEGER,ALLOCATABLE :: equiv(:,:,:), irr_k(:,:)
  REAL(DP),ALLOCATABLE :: eig1(:,:,:,:,:), eig2(:,:,:,:,:,:), &
  &                       nesting_irr(:,:), nesting_full(:,:,:,:)
  LOGICAL :: needwf = .FALSE.
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CALL mp_startup ()
  CALL environment_start ('NESTING')
  !
  ! ... Read pw.x input file and get prefix and outdir
  !
  CALL read_input_file ('PW', input_file_)
  !
  prefix_ = TRIM(prefix)
  tmp_dir = trimcheck(outdir)
  !
  ! ... Read XML file generated by pw.x
  !
  CALL read_file_new( needwf)
  !
  ! ... Number of k and spin for each magnetic treatment
  !
  IF (nspin == 2) THEN
     ns = 2
     IF(two_fermi_energies) THEN
        ef1 = ef_up
        ef2 = ef_dw
     ELSE
        ef1 = ef
        ef2 = ef
     END IF
  ELSE
     ns = 1
  END IF
  nk = nks / ns
  !
  ! ... Find equivalent k point in irr-BZ for whole BZ
  !
  ALLOCATE(equiv(nk1, nk2, nk3))
  CALL rotate_k_fs(equiv)
  !
  ! k in irreducible BZ
  !
  nirr_k = nk
  nks = nk1 * nk2 * nk3
  nkstot = nks
  ALLOCATE(nesting_irr(nirr_k,ns), nesting_full(nk1,nk2,nk3,ns), irr_k(3,nirr_k))
  !
  DO ik = 1, nirr_k
     irr_k(1:3,ik) = NINT(MATMUL(xk(1:3,ik), at(1:3,1:3)) * REAL((/nk1, nk2, nk3/), kind=dp))
     write(*,*) "debug ", irr_k(1:3,ik)
  END DO
  !
  ik = 0
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        DO i3 = 1, nk3
           ik = ik + 1
           xk(1:3,ik) = REAL((/i1, i2, i3/), DP) / REAL((/nk1, nk2, nk3/), DP)
           WHERE((/i1, i2, i3/)*2 >= (/nk1, nk2, nk3/)) xk(1:3,ik) = xk(1:3,ik) - 1.0_dp
           xk(1:3,ik) = MATMUL(bg(1:3,1:3), xk(1:3,ik))
           wk(ik) = 1.0_dp / REAL(nks, DP)
        END DO
     END DO
  END DO
  !
  ! Use tetrahedron method without symmetry
  !
  tetra_type = 2
  s_dummy(1:3,1:3,1:48) = 0
  t_rev_dummy(1:48) = 0
  DO i1 = 1, 3
     s_dummy(i1,i1,1:48) = 1
  END DO
  CALL opt_tetra_init(1, s_dummy, .False., t_rev_dummy, at, bg, npk, &
  &                   0, 0, 0, nk1, nk2, nk3, nks, xk, 1)
  !
  ! Nesting function function delta(e_k)*delta(e_{k+q})
  !
  nbnd_fs = b_high - b_low + 1
  ALLOCATE(eig1(nbnd_fs, nk3, nk2, nk1, ns), &
  &        eig2(nbnd_fs, nirr_k, nk3, nk2, nk1, ns))
  !
  ! ... Map e_k into whole BZ (Measured from E_F)
  !
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        DO i3 = 1, nk3
           IF(nspin == 2) THEN
              eig1(1:nbnd_fs,i3,i2,i1,1) = et(b_low:b_high, equiv(i1,i2,i3)     ) - ef1
              eig1(1:nbnd_fs,i3,i2,i1,2) = et(b_low:b_high, equiv(i1,i2,i3) + nk) - ef2
           ELSE
              eig1(1:nbnd_fs,i3,i2,i1,1) = et(b_low:b_high, equiv(i1,i2,i3)     ) - ef
           END IF
        END DO
     END DO
  END DO
  !
  ! Compute E(k+q) by shifting grid
  !
  DO ik = 1, nirr_k
     DO i1 = 1, nk1
        DO i2 = 1, nk2
           DO i3 = 1, nk3
              ikqv(1:3) = (/i1, i2, i3/) - 1 + irr_k(1:3,ik)
              ikqv(1:3) = MODULO(ikqv(1:3), (/nk1, nk2, nk3/)) + 1
              eig2(1:nbnd_fs,ik,i3,i2,i1,1:ns) = eig1(1:nbnd_fs,ikqv(3),ikqv(2),ikqv(1),1:ns)
           END DO
        END DO
     END DO
  END DO
  !
  nesting_irr(1:nirr_k,1:ns) = 0.0_dp
  !
  CALL nesting_delta1(nbnd_fs,ns,nirr_k,irr_k,eig1,eig2,nesting_irr)
  !
  ! Extend nesting factor to the full BZ
  !
  DO i3 = 1, nk3
     DO i2 = 1, nk2
        DO i1 = 1, nk1
           nesting_full(i1,i2,i3,1:ns) = nesting_irr(equiv(i1,i2,i3), 1:ns)
        END DO
     END DO
  END DO
  !
  ! ... Output in the FermiSurfer format
  !
  b_low = 1
  b_high = 1
  IF (nspin == 2) THEN
     !
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 1), TRIM(prefix) // "_nesting_delta1.frmsf")
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 2), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 2), "nesting_delta2.frmsf")
     !
  ELSE
     !
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 1), TRIM(prefix) // "_nesting_delta.frmsf")
     !
  END IF
  !
  ! Nesting function with (f_k - f_{k+q})/(e_{k+q} - e_k)
  !
  nesting_irr(1:nirr_k,1:ns) = 0.0_dp 
  CALL nesting_theta1(nbnd_fs,ns,nirr_k,eig1,eig2,nesting_irr)
  nesting_irr(1:nirr_k, 1:ns) = 2.0_dp * nesting_irr(1:nirr_k, 1:ns)
  IF (nspin == 1) nesting_irr(1:nirr_k, 1:ns) = 2.0_dp * nesting_irr(1:nirr_k, 1:ns)
  !
  ! Drude term
  !
  edos = 0.0_dp
  CALL opt_tetra_dos_t(eig1, nspin, nbnd_fs, nks, edos, dosf)
  DO ik = 1, nirr_k
    IF(ALL(irr_k(1:3,ik) == 0)) THEN
      nesting_irr(ik, 1:ns) = nesting_irr(ik, 1:ns) + dosf(1:ns)
    END IF
  END DO
  !
  ! Extend nesting factor to the full BZ
  !
  DO i3 = 1, nk3
     DO i2 = 1, nk2
        DO i1 = 1, nk1
           nesting_full(  i1,i2,i3,1:ns) = nesting_irr(equiv(i1,i2,i3), 1:ns)
        END DO
     END DO
  END DO
  !
  ! ... Output in the FermiSurfer format
  !
  b_low = 1
  b_high = 1
  IF (nspin == 2) THEN
     !
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      TRIM(tmp_dir) // TRIM(prefix) // "_nesting_chi1.frmsf")
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 2), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 2), &
     &                      TRIM(tmp_dir) // TRIM(prefix) // "_nesting_chi2.frmsf")
     !
  ELSE
     !
     CALL write_fermisurfer(nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      nesting_full(1:nk1, 1:nk2, 1:nk3, 1), &
     &                      TRIM(tmp_dir) // TRIM(prefix) // "_nesting_chi.frmsf")
     !
  END IF
  !
  DEALLOCATE(eig1, eig2, nesting_full, nesting_irr, equiv)
  !
  CALL environment_end ('NESTING')
#if defined(__MPI)
  CALL mp_global_end ( )
#endif
  !
END PROGRAM nesting
