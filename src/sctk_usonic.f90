!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!>
!> Module for the ultrasonic attenuation calculation.
!>
MODULE sctk_usonic
  !
  IMPLICIT NONE
  !
CONTAINS
!>
!> Compute the Fermi verocity
!>
SUBROUTINE calc_fvel()
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : alat, at
  USE fermisurfer_common, ONLY : b_low, b_high
  USE klist, ONLY : nks
  USE wvfct, ONLY : et
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE start_k, ONLY : nk1, nk2, nk3
  !
  USE sctk_val, ONLY : Fvel
  !
  IMPLICIT NONE
  !
  INTEGER :: ib, ik, ikp(3), ikm(3), nks0, nks1
  REAL(dp) :: de(3)
  !
  ALLOCATE(Fvel(b_low:b_high, nks, 3))
  Fvel(b_low:b_high,1:nks, 1:3) = 0.0_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(alat,nks0,nks1,nk3,nk2,nk1,b_low,b_high,at,Fvel,et,nks) &
  !$OMP & PRIVATE(ik,ib,ikp,ikm,de)
  !
  !$OMP DO
  DO ik = nks0, nks1
     DO ib = b_low, b_high
        !
        ikp(1) = MODULO(ik - 1 + nk2*nk3, nks) + 1
        ikp(2) = MODULO(ik - 1 + nk3,     nks) + 1
        ikp(3) = MODULO(ik - 1 + 1,       nks) + 1
        !
        ikm(1) = MODULO(ik - 1 - nk2*nk3, nks) + 1
        ikm(2) = MODULO(ik - 1 - nk3,     nks) + 1
        ikm(3) = MODULO(ik - 1 - 1,       nks) + 1
        !
        de(1:3) = et(ib, ikp(1:3)) - et(ib, ikm(1:3))
        de(1:3) = 0.5_dp * de(1:3) * REAL((/nk1,nk2,nk3/), dp)
        de(1:3) = MATMUL(at(1:3,1:3), de(1:3)) * alat
        Fvel(ib,ik,1) = de(1)
        Fvel(ib,ik,2) = de(1) + de(2)
        Fvel(ib,ik,3) = de(3)
        !
     END DO
     !
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  !
  CALL mp_sum( Fvel, world_comm )
  !
END SUBROUTINE calc_fvel
!>
!> Calc. ultrasonic attenuation with the tetrahedron method
!>
SUBROUTINE calc_usonic()
  !
  USE kinds, ONLY : DP
  USE ktetra, ONLY : wlsm, tetra, ntetra
  USE io_global, ONLY : stdout
  USE fermisurfer_common, ONLY : b_low, b_high
  USE klist, ONLY : nks
  USE wvfct, ONLY : et
  USE ener, ONLY : ef
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  !
  USE sctk_val, ONLY : dltF, Fvel
  !
  IMPLICIT NONE
  !
  INTEGER :: it, ib, ii, itetra(4), ntetra0, ntetra1
  REAL(dp) :: ei1(4,b_low:b_high), Fvel1(4,3,b_low:b_high), dl1(4,b_low:b_high), &
  &          Fvel2(3,3), dl2(3), &
  &          e(4), a(4,4), V, usonic0(2,3), usonic(2,3), tsmall(3,4)
  REAL(dp) :: dmin, dmax
  !
  WRITE(stdout,*) MINVAL(dltF(1,b_low:b_high,1:nks)), MAXVAL(dltF(1,b_low:b_high,1:nks)), &
  &                SUM(dltF(1,b_low:b_high,1:nks)) / REAL((b_high-b_low+1) * nks, dp)
  usonic(1:2,1:3) = 0.0_dp
  !
  CALL divide(world_comm, ntetra,ntetra0,ntetra1)
  !
  dmax = -1e10_dp
  dmin = 1e10_dp
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(b_low,b_high,wlsm,Fvel,dltF,usonic,dmin,dmax,et,tetra,ef,ntetra0,ntetra1) &
  !$OMP & PRIVATE(it,ib,ii,a,e,ei1,Fvel1,dl1,Fvel2,dl2,V,usonic0,itetra,tsmall)
  !
  !$OMP DO REDUCTION(+:usonic) REDUCTION(max:dmax) REDUCTION(min:dmin)
  DO it = ntetra0, ntetra1
     !
     ei1(  1:4,      b_low:b_high) = 0.0_dp
     Fvel1(1:4, 1:3, b_low:b_high) = 0.0_dp
     dl1(  1:4,      b_low:b_high) = 0.0_dp
     !
     DO ii = 1, 20
        !
        DO ib = b_low, b_high
           ei1(  1:4,    ib) = MATMUL(wlsm(1:4,1:20), et(     ib, tetra(1:20, it))) - ef
           Fvel1(1:4,1:3,ib) = MATMUL(wlsm(1:4,1:20), Fvel(   ib, tetra(1:20, it), 1:3))
           dl1(  1:4,    ib) = MATMUL(wlsm(1:4,1:20), dltF(1, ib, tetra(1:20, it)))
        END DO
     END DO
     !
     DO ib = b_low, b_high
        !
        itetra(1) = 0
        e(1:4) = ei1(1:4,ib)
        CALL hpsort (4, e, itetra)
        !
        DO ii = 1, 4
           a(ii,1:4) = (0.0_dp - e(1:4)) / (e(ii) - e(1:4))
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
           Fvel2(1:3,1:3) = MATMUL(tsmall(1:3,1:4), Fvel1(itetra(1:4), 1:3, ib))
           dl2(1:3) = MATMUL(tsmall(1:3,1:4), dl1(itetra(1:4), ib))
           !
           dmin = MIN(dmin, MINVAL(dl2(1:3)))
           dmax = MAX(dmax, MAXVAL(dl2(1:3)))
           ! 
           CALL usonic_step2(Fvel2,dl2,usonic0)
           usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
           !
        ELSE IF(e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) THEN
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
           Fvel2(1:3,1:3) = MATMUL(tsmall(1:3,1:4), Fvel1(itetra(1:4), 1:3, ib))
           dl2(1:3) = MATMUL(tsmall(1:3,1:4), dl1(itetra(1:4), ib))
           ! 
           dmin = MIN(dmin, MINVAL(dl2(1:3)))
           dmax = MAX(dmax, MAXVAL(dl2(1:3)))
           !
           CALL usonic_step2(Fvel2,dl2,usonic0)
           usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
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
           Fvel2(1:3,1:3) = MATMUL(tsmall(1:3,1:4), Fvel1(itetra(1:4), 1:3, ib))
           dl2(1:3) = MATMUL(tsmall(1:3,1:4), dl1(itetra(1:4), ib))
           ! 
           dmin = MIN(dmin, MINVAL(dl2(1:3)))
           dmax = MAX(dmax, MAXVAL(dl2(1:3)))
           !
           CALL usonic_step2(Fvel2,dl2,usonic0)
           usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
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
           Fvel2(1:3,1:3) = MATMUL(tsmall(1:3,1:4), Fvel1(itetra(1:4), 1:3, ib))
           dl2(1:3) = MATMUL(tsmall(1:3,1:4), dl1(itetra(1:4), ib))
           !
           dmin = MIN(dmin, MINVAL(dl2(1:3)))
           dmax = MAX(dmax, MAXVAL(dl2(1:3)))
           !
           CALL usonic_step2(Fvel2,dl2,usonic0)
           usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
           !
        END IF
        !
     END DO ! ib
     !
  END DO ! it
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  CALL mp_sum( usonic, world_comm )
  !
  WRITE(stdout,*) "Min & Max : ", dmin, dmax
  WRITE(stdout,*) ""
  WRITE(stdout,*) "Ultrasonic attenuation coefficient"
  WRITE(stdout,*) "q, a_N, a_S, a_S / a_N"
  WRITE(stdout,*) "100 : ", usonic(1,1), usonic(2,1), usonic(2,1) / usonic(1,1)
  WRITE(stdout,*) "110 : ", usonic(1,2), usonic(2,2), usonic(2,2) / usonic(1,2)
  WRITE(stdout,*) "001 : ", usonic(1,3), usonic(2,3), usonic(2,3) / usonic(1,3)
  WRITE(stdout,*) ""
  !
END SUBROUTINE calc_usonic
!>
!> 2nd step of tetrahedra method.
!>
SUBROUTINE usonic_step2(Fvel1,dl1,usonic)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: Fvel1(3,3) !< The Fermi velocity
  REAL(dp),INTENT(IN) :: dl1(3) !< @f$\Delta_{n k}@f$
  REAL(dp),INTENT(OUT) :: usonic(2,3) !< ultrasonic attenuation
  !
  INTEGER :: id, ii, itri(3)
  REAL(dp) :: e(3), a(3,3), dl2(2), V, usonic0
  !
  usonic(1:2,1:3) = 0.0_dp
  !
  DO id = 1, 3
     !
     IF(MAXVAL(ABS(Fvel1(id,1:3))) < 1e-10_dp) THEN
        !$OMP CRITICAL   
        WRITE(*,*) Fvel1(id,1:3)
        stop "Fvel = 0 !!"
        !$OMP END CRITICAL
     END IF
     !
     itri(1) = 0
     e(1:3) = Fvel1(1:3, id)
     CALL hpsort (3, e, itri)
     !
     DO ii = 1, 3
        a(ii,1:3) = (0.0_dp - e(1:3)) / (e(ii) - e(1:3))
     END DO
     !
     IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
        !
        !V = 2.0_dp * a(2,1) * a(3,1) / (0.0_dp - e(1)) 
        V = 2.0_dp * a(2,1)           / (e(3) - e(1)) 
        !
        dl2(1) = dl1(itri(1)) * a(1,2) + dl1(itri(2)) * a(2,1)
        dl2(2) = dl1(itri(1)) * a(1,3) + dl1(itri(3)) * a(3,1)
        !
        IF(ABS(dl2(1) - dl2(2)) < 1e-10_dp) THEN
           usonic0 = 0.5_dp * (1.0_dp - tanh(0.25_dp * ABS(sum(dl2(1:2)))))
        ELSE IF(dl2(1) * dl2(2) >= 0.0_dp) THEN
           usonic0 = (1.0_dp + tanh(0.5_dp * ABS(dl2(1)))) / (1.0_dp + tanh(0.5_dp * ABS(dl2(2))))
           usonic0 = LOG(usonic0) / (ABS(dl2(1)) - ABS(dl2(2)))
        ELSE
           usonic0 = (1.0_dp + tanh(0.5_dp * ABS(dl2(1)))) * (1.0_dp + tanh(0.5_dp * ABS(dl2(2))))
           usonic0 = LOG(usonic0) / (ABS(dl2(1)) + ABS(dl2(2)))
        END IF
        !
        usonic(1,id) = V
        usonic(2,id) = V * 2.0_dp * usonic0
        !
     ELSE IF((e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) .OR. (e(2) <= 0.0_dp .AND. 0.0_dp < e(3))) THEN
        !
        !V = 2.0_dp * a(1,3) * a(2,3) / (e(3) - 0.0_dp)
        V = 2.0_dp * a(1,3)           / (e(3) - e(2))
        !
        dl2(1) = dl1(itri(3)) * a(3,1) + dl1(itri(1)) * a(1,3)
        dl2(2) = dl1(itri(3)) * a(3,2) + dl1(itri(2)) * a(2,3)
        !
        IF(ABS(dl2(1) - dl2(2)) < 1e-10_dp) THEN
           usonic0 = 0.5_dp * (1.0_dp - tanh(0.25_dp * ABS(sum(dl2(1:2)))))
        ELSE IF(dl2(1) * dl2(2) >= 0.0_dp) THEN
           usonic0 = (1.0_dp + tanh(0.5_dp * ABS(dl2(1)))) / (1.0_dp + tanh(0.5_dp * ABS(dl2(2))))
           usonic0 = LOG(usonic0) / (ABS(dl2(1)) - ABS(dl2(2)))
        ELSE
           usonic0 = (1.0_dp + tanh(0.5_dp * ABS(dl2(1)))) * (1.0_dp + tanh(0.5_dp * ABS(dl2(2))))
           usonic0 = LOG(usonic0) / (ABS(dl2(1)) + ABS(dl2(2)))
        END IF
        !
        usonic(1,id) = V
        usonic(2,id) = V * 2.0_dp * usonic0
        !
     END IF
     !
  END DO ! id
  !
END SUBROUTINE usonic_step2
  !
END MODULE sctk_usonic
