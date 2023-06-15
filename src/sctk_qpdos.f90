!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_qpdos
  !
  IMPLICIT NONE
  !
CONTAINS
!>
!> Construct the energy grid for qpdos
!>
SUBROUTINE egrid()
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : RYTOEV
  !
  USE sctk_val, ONLY : e0, emax, ne, nx, sdos, xi0
  !
  IMPLICIT NONE
  !
  INTEGER :: ie
  REAL(dp) :: de
  !
  emax = emax / (RYTOEV * 1.0e3_dp)
  nx = ne * 2 - 1
  !
  DEALLOCATE(xi0)
  ALLOCATE(e0(ne), sdos(ne), xi0(nx))
  de = emax / REAL(ne, dp)
  !
  DO ie = 1, ne
     e0(ie) = de * REAL(ie, dp)
  END DO
  !
  DO ie = 1, nx
     xi0(ie) = de * REAL(ie - ne, dp)
  END DO
  !
  sdos(1:ne) = 0.0_dp
  !
END SUBROUTINE egrid
!>
!> Calc. qpdos by using 4D interpolated integration method
!>
SUBROUTINE calc_sdos()
  !
  USE kinds, ONLY : DP
  USE start_k, ONLY : nk1, nk2, nk3
  USE klist, ONLY : nks
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE fermisurfer_common, ONLY : b_low
  !
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE fermisurfer_common,   ONLY : b_low, b_high
  USE sctk_val, ONLY : dltF, e0, ne, nx, sdos
  USE sctk_tetra, ONLY : calc_dosk
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, ix, ie, ip, ii, ik2, ix2, ikv1(3), iv4d(4,5,24), ipenta(5), nks0, nks1
  REAL(dp) :: e(5), dos2(5), a(5,5), dosk(nx,elph_nbnd_min:elph_nbnd_max,nks)
  !
  CALL calc_dosk(dosk)
  CALL mp_sum(dosk, world_comm )
  !
  CALL penta_4d(iv4d)
  CALL divide(world_comm, nks,nks0,nks1)
  !
  sdos = 0.0_dp
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,b_low,b_high,nx,ne,nk1,nk2,nk3,iv4d,dltF,dosk,e0,sdos) &
  !$OMP & PRIVATE(ik,ib,ix,ip,ii,ie,ik2,ix2,e,a,dos2,ikv1,ipenta)
  !
  !$OMP DO REDUCTION(+: sdos)
  DO ik = nks0, nks1
     DO ib = b_low, b_high
        DO ix = 1, nx - 1
           !
           DO ip = 1, 24
              !
              DO ii = 1, 5
                 !
                 ikv1(1) = (ik - 1) / (nk2*nk3)
                 ikv1(2) = (ik - 1 - ikv1(1)*nk2*nk3) / nk3
                 ikv1(3) =  ik - 1 - ikv1(1)*nk2*nk3 - ikv1(2)*nk3
                 ikv1(1:3) = ikv1(1:3) + iv4d(1:3,ii,ip)
                 ikv1(1:3) = MODULO(ikv1(1:3), (/nk1, nk2, nk3/))
                 ik2 = 1 + ikv1(3) + nk3*ikv1(2) + nk3*nk2*ikv1(1)
                 !
                 ix2 = ix + iv4d(4,ii,ip)
                 !
                 e(   ii) = dltF(ix2,ib,ik2)
                 dos2(ii) = dosk(ix2,ib,ik2)
                 !
              END DO
              !
              ipenta(1) = 0
              CALL hpsort (5, e, ipenta)
              !
              DO ie = 1, ne
                 !
                 DO ii = 1, 5
                    a(1:5,ii) = (e0(ie) - e(ii)) / (e(1:5) - e(ii))
                 END DO
                 !
                 IF(e(1) < e0(ie) .AND. e0(ie) <= e(2)) THEN
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(2,1) * a(3,1) * a(4,1) * a(5,1) / (e0(ie) - e(1)) &
                    &        * ( dos2(ipenta(1)) * a(1,2) + dos2(ipenta(2)) * a(2,1) & 
                    &          + dos2(ipenta(1)) * a(1,3) + dos2(ipenta(3)) * a(3,1) &
                    &          + dos2(ipenta(1)) * a(1,4) + dos2(ipenta(4)) * a(4,1) &
                    &          + dos2(ipenta(1)) * a(1,5) + dos2(ipenta(5)) * a(5,1) )
                    !
                 ELSE IF(e(2) < e0(ie) .AND. e0(ie) <= e(3)) THEN
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(3,1) * a(4,1) * a(5,1) * a(2,3) / (e0(ie) - e(1)) &
                    &        * ( dos2(ipenta(1)) * a(1,3) + dos2(ipenta(3)) * a(3,1) &
                    &          + dos2(ipenta(1)) * a(1,4) + dos2(ipenta(4)) * a(4,1) &
                    &          + dos2(ipenta(1)) * a(1,5) + dos2(ipenta(5)) * a(5,1) &
                    &          + dos2(ipenta(2)) * a(2,3) + dos2(ipenta(3)) * a(3,2) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(4,1) * a(5,1) * a(3,2) * a(2,4) / (e0(ie) - e(1)) &
                    &        * ( dos2(ipenta(1)) * a(1,4) + dos2(ipenta(4)) * a(4,1) &
                    &          + dos2(ipenta(1)) * a(1,5) + dos2(ipenta(5)) * a(5,1) &
                    &          + dos2(ipenta(2)) * a(2,3) + dos2(ipenta(3)) * a(3,2) &
                    &          + dos2(ipenta(2)) * a(2,4) + dos2(ipenta(4)) * a(4,2) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(5,1) * a(3,2) * a(4,2) * a(2,5) / (e0(ie) - e(1)) &
                    &        * ( dos2(ipenta(1)) * a(1,5) + dos2(ipenta(5)) * a(5,1) &
                    &          + dos2(ipenta(2)) * a(2,3) + dos2(ipenta(3)) * a(3,2) &
                    &          + dos2(ipenta(2)) * a(2,4) + dos2(ipenta(4)) * a(4,2) & 
                    &          + dos2(ipenta(2)) * a(2,5) + dos2(ipenta(5)) * a(5,2) )
                    !
                ELSE IF(e(3) < e0(ie) .AND. e0(ie) <= e(4)) THEN
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(3,5) * a(2,5) * a(1,5) * a(4,3) / (e(5) - e0(ie)) &
                    &        * ( dos2(ipenta(5)) * a(5,3) + dos2(ipenta(3)) * a(3,5) &
                    &          + dos2(ipenta(5)) * a(5,2) + dos2(ipenta(2)) * a(2,5) &
                    &          + dos2(ipenta(5)) * a(5,1) + dos2(ipenta(1)) * a(1,5) &
                    &          + dos2(ipenta(4)) * a(4,3) + dos2(ipenta(3)) * a(3,4) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(2,5) * a(1,5) * a(3,4) * a(4,2) / (e(5) - e0(ie)) &
                    &        * ( dos2(ipenta(5)) * a(5,2) + dos2(ipenta(2)) * a(2,5) &
                    &          + dos2(ipenta(5)) * a(5,1) + dos2(ipenta(1)) * a(1,5) &
                    &          + dos2(ipenta(4)) * a(4,3) + dos2(ipenta(3)) * a(3,4) &
                    &          + dos2(ipenta(4)) * a(4,2) + dos2(ipenta(2)) * a(2,4) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(1,5) * a(3,4) * a(2,4) * a(4,1) / (e(5) - e0(ie)) &
                    &        * ( dos2(ipenta(5)) * a(5,1) + dos2(ipenta(1)) * a(1,5) &
                    &          + dos2(ipenta(4)) * a(4,3) + dos2(ipenta(3)) * a(3,4) &
                    &          + dos2(ipenta(4)) * a(4,2) + dos2(ipenta(2)) * a(2,4) &
                    &          + dos2(ipenta(4)) * a(4,1) + dos2(ipenta(1)) * a(1,4) )
                    !
                 ELSE IF(e(4) < e0(ie) .AND. e0(ie) <= e(5)) THEN
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(4,5) * a(3,5) * a(2,5) * a(1,5) / (e(5) - e0(ie)) &
                    &        * ( dos2(ipenta(5)) * a(5,4) + dos2(ipenta(4)) * a(4,5) &
                    &          + dos2(ipenta(5)) * a(5,3) + dos2(ipenta(3)) * a(3,5) &
                    &          + dos2(ipenta(5)) * a(5,2) + dos2(ipenta(2)) * a(2,5) &
                    &          + dos2(ipenta(5)) * a(5,1) + dos2(ipenta(1)) * a(1,5) )
                    !
                 END IF
                 !
              END DO ! ie
              !
           END DO ! ip
           !
        END DO ! ix
     END DO ! ib
  END DO ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  sdos(1:ne) = sdos(1:ne) / 24.0_dp * (e0(2) - e0(1))
  !
  CALL mp_sum( sdos, world_comm )
  !
END SUBROUTINE calc_sdos
!>
!> Define corners of the hyper pentahedron
!>
SUBROUTINE penta_4d(iv4d)
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(OUT) :: iv4d(4,5,24)
  !
  INTEGER :: i1, i2, i3, i4, ip
  !
  ip = 0
  DO i1 = 1, 4
     DO i2 = 1, 4
        IF(i2 == i1) CYCLE
        DO i3 = 1, 4
           IF(i3 == i1 .OR. i3 == i2) CYCLE
           DO i4 = 1, 4
              IF(i4 == i1 .OR. i4 == i2 .OR. i4 == i3) CYCLE
              !
              ip = ip + 1
              !
              iv4d(1:4,1,ip) = 0
              !
              iv4d(1:4,2,ip) = iv4d(1:4,1,ip)
              iv4d( i1,2,ip) = iv4d( i1,2,ip) + 1
              !
              iv4d(1:4,3,ip) = iv4d(1:4,2,ip)
              iv4d( i2,3,ip) = iv4d( i2,3,ip) + 1
              !
              iv4d(1:4,4,ip) = iv4d(1:4,3,ip)
              iv4d( i3,4,ip) = iv4d( i3,4,ip) + 1
              !
              iv4d(1:4,5,ip) = iv4d(1:4,4,ip)
              iv4d( i4,5,ip) = iv4d( i4,5,ip) + 1
              !
           END DO
        END DO
     END DO
  END DO
  !
END SUBROUTINE penta_4d
  !
END MODULE sctk_qpdos
