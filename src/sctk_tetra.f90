!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_tetra
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Bi-linear interpolation
!
SUBROUTINE interpol_indx(ng,kvec,kintp,wintp)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: ng(3)
  REAL(DP),INTENT(in) :: kvec(3)
  INTEGER,INTENT(out) :: kintp(8)
  REAL(DP),INTENT(out) :: wintp(8)
  !
  INTEGER :: ii
  REAL(DP) :: xv(3)
  integer :: ikv0(3), ikv1(3), i1, i2, i3
  !
  xv(1:3) = kvec(1:3) * real(ng(1:3), dp)
  ikv0(1:3) = floor(xv(1:3))
  xv(1:3) = xv(1:3) - real(ikv0(1:3), dp)
  !
  ii = 0
  do i1 = 0, 1
     do i2 = 0, 1
        do i3 = 0, 1
           !
           ii = ii + 1
           !
           ikv1(1:3) = ikv0(1:3) + (/i1, i2, i3/)
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           kintp(ii) = 1 + ikv1(3) + ng(3) * ikv1(2) + ng(3) * ng(2) * ikv1(1)
           !
           wintp(ii) = xv(1)**i1 * (1.0_dp - xv(1))**(1 - i1) &
           &         * xv(2)**i2 * (1.0_dp - xv(2))**(1 - i2) &
           &         * xv(3)**i3 * (1.0_dp - xv(3))**(1 - i3) 
           !
        end do
     end do
  end do
  !
END SUBROUTINE interpol_indx
!
! Compute DOS at each k flagment
!
SUBROUTINE calc_dosk(dosd)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE klist, ONLY : nks
  USE wvfct, ONLY : et
  USE ktetra, ONLY : wlsm, tetra, ntetra
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE ener, ONLY : ef
  !
  USE sctk_val, ONLY : nx, xi0
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(OUT) :: dosd(nx,elph_nbnd_min:elph_nbnd_max,nks)
  !
  INTEGER :: ibnd, it, ii, ix, itetra(4), ntetra0, ntetra1
  REAL(dp) :: ei(4,elph_nbnd_min:elph_nbnd_max), e(4), a(4,4), w1(nx,4), V
  !
  dosd(1:nx, elph_nbnd_min:elph_nbnd_max, 1:nks) = 0.0_dp
  !
  CALL divide(world_comm, ntetra, ntetra0, ntetra1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ntetra0,ntetra1,elph_nbnd_min,elph_nbnd_max,nx,wlsm,et,ef,xi0,dosd,tetra) &
  !$OMP & PRIVATE(it,ii,ibnd,ix,ei,w1,a,e,V,itetra)
  !
  DO it = ntetra0, ntetra1
     !
     DO ibnd = elph_nbnd_min, elph_nbnd_max
        ei(1:4, ibnd) = MATMUL(wlsm(1:4,1:20), et(ibnd,tetra(1:20, it))) - ef
     END DO
     !
     !$OMP DO
     DO ibnd = elph_nbnd_min, elph_nbnd_max
        !
        itetra(1) = 0
        e(1:4) = ei(1:4,ibnd)
        CALL hpsort (4, e, itetra)
        !
        DO ix = 1, nx
           !
           DO ii = 1, 4
              a(ii,1:4) = (xi0(ix) - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) < xi0(ix) .AND. xi0(ix) <= e(2)) THEN
              !
              V = a(2,1) * a(3,1) * a(4,1) / (xi0(ix) - e(1))
              !
              w1(ix,itetra(1)) = a(1,2) + a(1,3) + a(1,4)
              w1(ix,itetra(2:4)) = a(2:4,1)
              w1(ix,1:4) = w1(ix,1:4) * V
              !
           ELSE IF(e(2) < xi0(ix) .AND. xi0(ix) <= e(3)) THEN
              !
              V = a(2,3) * a(3,1) + a(3,2) * a(2,4)
              !
              w1(ix,itetra(1)) = a(1,4) * V + a(1,3) * a(3,1) * a(2,3)
              w1(ix,itetra(2)) = a(2,3) * V + a(2,4) * a(2,4) * a(3,2)
              w1(ix,itetra(3)) = a(3,2) * V + a(3,1) * a(3,1) * a(2,3)
              w1(ix,itetra(4)) = a(4,1) * V + a(4,2) * a(2,4) * a(3,2)
              !
              V = 1.0_dp / (e(4) - e(1))
              w1(ix,1:4) = w1(ix,1:4) * V
              !
           ELSE IF(e(3) < xi0(ix) .AND. xi0(ix) < e(4)) THEN
              !
              V = a(1,4) * a(2,4) * a(3,4) / (e(4) - xi0(ix))
              w1(ix,itetra(1:3)) = a(1:3,4)
              w1(ix,itetra(4)) = a(4,1) + a(4,2) + a(4,3)
              w1(ix,1:4) = w1(ix,1:4) * V
              !
           ELSE
              !
              w1(ix,1:4) = 0.0_dp
              !
           END IF
           !
        END DO ! ix
        !
        dosd(1:nx, ibnd, tetra(1:20, it)) = dosd(1:nx,ibnd,    tetra(1:20, it)) &
        &                            + MATMUL(w1(1:nx,1:4), wlsm(1:4,1:20))
        !
     END DO ! ibnd
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !$OMP END PARALLEL
  !
  dosd(  1:nx, elph_nbnd_min:elph_nbnd_max, 1:nks) = &
  & dosd(1:nx, elph_nbnd_min:elph_nbnd_max, 1:nks) / REAL(ntetra, dp)
  !
END SUBROUTINE calc_dosk
!
! Integration weight with tetrahedron method
!
SUBROUTINE tetraweight(wghtd)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE wvfct, ONLY : et
  USE ktetra, ONLY : wlsm, tetra, ntetra
  USE klist, ONLY : nks
  USE ener, ONLY : ef
  USE lsda_mod, ONLY : nspin
  !
  USE sctk_val, ONLY : nmf, nb, bdsp
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(dp),INTENT(OUT) :: wghtd((nmf+1)*nb(2),nb(1),nks)
  !
  INTEGER :: it, ibnd, ii, itetra(4)
  REAL(dp) :: thr = 1e-8_dp, V
  REAL(dp) :: e(4), a(4,4), ei0(4,nb(1)), ej0(4,nb(2)), ei1(4), ej1(4,nb(2)), tsmall(4,4)
  COMPLEX(dp) :: w1(nb(2)*(nmf+1),4), w2(nb(2)*(nmf+1),4)
  !
  wghtd(1:nb(2)*(nmf+1),1:nb(1),1:nks) = 0.0_dp
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(wlsm,nmf,et,wghtd,thr,ef,tetra,nks,nb,bdsp,ntetra) &
  !$OMP & PRIVATE(it,ii,ibnd,itetra,ei0,ej0,ei1,ej1,w1,w2,a,e,V,tsmall)
  !
  DO it = 1, ntetra
     !
     DO ibnd = 1, nb(1)
        ei0(1:4,ibnd) = MATMUL(wlsm(1:4,1:20), et(ibnd+bdsp(1), tetra(1:20,it)      )) - ef
     END DO
     DO ibnd = 1, nb(2)
        ej0(1:4,ibnd) = MATMUL(wlsm(1:4,1:20), et(ibnd+bdsp(2), tetra(1:20,it) + nks)) - ef
     END DO
     !
     !$OMP DO
     DO ibnd = 1, nb(1)
        !
        w1(1:(nmf+1)*nb(2),1:4) = 0.0_dp
        !
        itetra(1) = 0
        e(1:4) = ei0(1:4, ibnd)
        CALL hpsort (4, e, itetra)
        !
        DO ii = 1, 4
           a(ii,1:4) = (0.0_dp - e(1:4)) / (e(ii) - e(1:4))
        END DO
        !
        IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2) ) THEN
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
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
              ei1(1:4        ) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4), ibnd))
              ej1(1:4,1:nb(2)) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nb(2)))
              !
              CALL tetra2(ei1,ej1,w2)
              !
              w1(1:(nmf+1)*nb(2),itetra(1:4)) = w1(1:(nmf+1)*nb(2),          itetra(1:4)) &
              &                    + V * MATMUL(w2(1:(nmf+1)*nb(2),1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0.0_dp ) THEN
           !
           ! D - 1
           !
           ei1(1:4        ) = ei0(1:4, ibnd)
           ej1(1:4,1:nb(2)) = ej0(1:4, 1:nb(2))
           !
           CALL tetra2(ei1,ej1,w2)
           !
           w1(1:(nmf+1)*nb(2),1:4) = w1(1:(nmf+1)*nb(2),1:4) &
           &                       + w2(1:(nmf+1)*nb(2),1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        wghtd(1:(nmf+1)*nb(2),ibnd,tetra(1:20,it)) = wghtd(1:(nmf+1)*nb(2),ibnd,  tetra(1:20,it)) &
        &                                      + MATMUL(w1(1:(nmf+1)*nb(2),1:4), wlsm(1:4,1:20))
        !
     END DO ! ibnd
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  IF(nspin == 1) THEN
     wghtd(1:(nmf+1)*nb(2),1:nb(1),1:nks) = wghtd(1:(nmf+1)*nb(2),1:nb(1),1:nks) &
     &                            * 4.0_dp / (REAL(ntetra, dp) * omega)
  ELSE
     wghtd(1:(nmf+1)*nb(2),1:nb(1),1:nks) = wghtd(1:(nmf+1)*nb(2),1:nb(1),1:nks) &
     &                            * 2.0_dp / (REAL(ntetra, dp) * omega)
  END IF
  !
END SUBROUTINE tetraweight
!
!-----------------------------------------------------------------------
SUBROUTINE tetra2(ei0,ej0,w0)
  !---------------------------------------------------------------------
  !
  ! This routine take the unoccupied region.
  !
  USE kinds, ONLY : dp
  !
  USE sctk_val, ONLY : nmf, nb
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei0(4), ej0(4,nb(2))
  COMPLEX(dp),INTENT(OUT) :: w0(0:nmf,nb(2),4)
  !
  INTEGER :: ii, ibnd, itetra(4)
  REAL(dp) :: V, de(4), thr = 1.0e-8_dp
  REAL(dp) :: e(4), a(4,4), tsmall(4,4)
  COMPLEX(dp) :: w1(0:nmf,4)
  !
  w0(0:nmf,1:nb(2),1:4) = 0.0_dp
  !
  DO ibnd = 1, nb(2)
     !
     e(1:4) = ej0(1:4,ibnd)
     itetra(1) = 0
     call hpsort (4, e, itetra)
     !
     DO ii = 1, 4
        a(ii,1:4) = ( 0.0_dp - e(1:4) ) / (e(ii) - e(1:4))
     END DO
     !
     IF(0_dp <= e(1)) THEN
        !
        ! A - 1
        !
        de(1:4) = e(1:4) - ei0(itetra(1:4))
        CALL lindhard(de,w1)
        w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) + w1(0:nmf,1:4)
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
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
           de(1:4) = MATMUL(tsmall(1:4,1:4), (e(1:4) - ei0(itetra(1:4))))
           CALL lindhard(de,w1)
           w0(0:nmf,ibnd,itetra(1:4)) = w0(0:nmf,ibnd,itetra(1:4)) &
           &               + V * MATMUL(w1(0:nmf,1:4), tsmall(1:4,1:4))
           !        
        END IF
        !
     END IF
     !
  END DO
  !
END SUBROUTINE tetra2
!
! Tetarahedra method for lindhard function
!
SUBROUTINE lindhard(de,w)
  !
  USE kinds, ONLY : DP
  USE dfpt_tetra_mod, ONLY : dfpt_tetra_lindhard_1234, dfpt_tetra_lindhard_1231, &
  &                          dfpt_tetra_lindhard_1233, dfpt_tetra_lindhard_1221, &
  &                          dfpt_tetra_lindhard_1222, dfpt_tetra_lindhard_1211
  !
  USE sctk_val, ONLY : nmf, mf
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: de(4)
  COMPLEX(dp),INTENT(OUT) :: w(0:nmf,4)
  !
  INTEGER :: ii, imf, itetra(4)
  REAL(dp) :: e(4), le(4), thr, thr2, em(4), w2(2,4)
  !
  ! Static part
  !
  itetra(1) = 0
  e(1:4) = de(1:4)
  call hpsort (4, e, itetra)
  !
  thr = MAXVAL(e(1:4)) * 1e-3_dp
  thr2 = 1.0e-12_dp
  !
  DO ii = 1, 4
     IF(e(ii) < thr2) THEN
        IF(ii == 3) THEN
           write(*,*) e, de, thr2
           CALL errore("lindhard", "Nesting occurs.", 1)
        END IF
        le(ii) = 0.0_dp
        e(ii) = 0.0_dp
     ELSE
        le(ii) = LOG(e(ii))
     END IF
  END DO
  !
  IF(ABS(e(4) - e(3)) < thr ) THEN
     IF(ABS(e(4) - e(2)) < thr ) THEN
        IF(ABS(e(4) - e(1)) < thr ) THEN
           !
           ! e(4) = e(3) = e(2) = e(1)
           !
           w(0,itetra(4)) = 0.25_dp / e(4)
           w(0,itetra(3)) = w(0,itetra(4))
           w(0,itetra(2)) = w(0,itetra(4))
           w(0,itetra(1)) = w(0,itetra(4))
           !
        ELSE
           !
           ! e(4) = e(3) = e(2)
           !
           w(0,itetra(4)) = dfpt_tetra_lindhard_1211(e(4),e(1),le(4),le(1))
           w(0,itetra(3)) = w(0,itetra(4))
           w(0,itetra(2)) = w(0,itetra(4))
           w(0,itetra(1)) = dfpt_tetra_lindhard_1222(e(1),e(4),le(1),le(4))
           !
           IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
              CALL errore("lindhard", "4=3=2", 1)
           END IF
           !
        END IF
     ELSE IF(ABS(e(2) - e(1)) < thr ) THEN
        !
        ! e(4) = e(3), e(2) = e(1)
        !
        w(0,itetra(4)) = dfpt_tetra_lindhard_1221(e(4),e(2), le(4),le(2))
        w(0,itetra(3)) = w(0,itetra(4))
        w(0,itetra(2)) = dfpt_tetra_lindhard_1221(e(2),e(4), le(2),le(4))
        w(0,itetra(1)) = w(0,itetra(2))
        !
        IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
           CALL errore("lindhard", "4=3 2=1", 1)
        END IF
        !
     ELSE
        !
        ! e(4) = e(3)
        !
        w(0,itetra(4)) = dfpt_tetra_lindhard_1231(e(4),e(1),e(2),le(4),le(1),le(2))
        w(0,itetra(3)) = w(0,itetra(4))
        w(0,itetra(2)) = dfpt_tetra_lindhard_1233(e(2),e(1),e(4),le(2),le(1),le(4))
        w(0,itetra(1)) = dfpt_tetra_lindhard_1233(e(1),e(2),e(4),le(1),le(2),le(4))
        !
        IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
           CALL errore("lindhard", "4=3", 1)
        END IF
        !
     END IF
  ELSE IF(ABS(e(3) - e(2)) < thr) THEN
     IF(ABS(e(3) - e(1)) < thr) THEN
        !
        ! e(3) = e(2) = e(1)
        !
        w(0,itetra(4)) = dfpt_tetra_lindhard_1222(e(4),e(3), le(4),le(3))
        w(0,itetra(3)) = dfpt_tetra_lindhard_1211(e(3),e(4), le(3),le(4))
        w(0,itetra(2)) = w(0,itetra(3))
        w(0,itetra(1)) = w(0,itetra(3))
        !
        IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
           CALL errore("lindhard", "3=2=1", 1)
        END IF
        !
     ELSE
        !
        ! e(3) = e(2)
        !
        w(0,itetra(4)) = dfpt_tetra_lindhard_1233(e(4),e(1),e(3),le(4),le(1),le(3))
        w(0,itetra(3)) = dfpt_tetra_lindhard_1231(e(3),e(1),e(4),le(3),le(1),le(4))
        w(0,itetra(2)) = w(0,itetra(3))
        w(0,itetra(1)) = dfpt_tetra_lindhard_1233(e(1),e(4),e(3),le(1),le(4),le(3))
        !
        IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
           CALL errore("lindhard", "3=2", 1)
        END IF
        !
     END IF
  ELSE IF(ABS(e(2) - e(1)) < thr) THEN
     !
     ! e(2) = e(1)
     !
     w(0,itetra(4)) = dfpt_tetra_lindhard_1233(e(4),e(3),e(2),le(4),le(3),le(2))
     w(0,itetra(3)) = dfpt_tetra_lindhard_1233(e(3),e(4),e(2),le(3),le(4),le(2))
     w(0,itetra(2)) = dfpt_tetra_lindhard_1231(e(2),e(3),e(4),le(2),le(3),le(4))
     w(0,itetra(1)) = w(0,itetra(2))
     !
     IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
        CALL errore("lindhard", "2=1", 1)
     END IF
     !
  ELSE
     !
     ! DIFferent each other.
     !
     w(0,itetra(4)) = dfpt_tetra_lindhard_1234(e(4),e(1),e(2),e(3),le(4),le(1),le(2),le(3))
     w(0,itetra(3)) = dfpt_tetra_lindhard_1234(e(3),e(1),e(2),e(4),le(3),le(1),le(2),le(4))
     w(0,itetra(2)) = dfpt_tetra_lindhard_1234(e(2),e(1),e(3),e(4),le(2),le(1),le(3),le(4))
     w(0,itetra(1)) = dfpt_tetra_lindhard_1234(e(1),e(2),e(3),e(4),le(1),le(2),le(3),le(4))
     !
     IF(ANY(REAL(w(0,itetra(1:4)), dp) < 0.0_dp)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') REAL(w(0,itetra(1:4)), dp)
        CALL errore("lindhard", "Something wrong.", 1)
     END IF
     !
  END IF
  !
  ! w /= 0 part
  !
  DO imf = 1, nmf
     !
     em(1:4) = e(1:4) / mf(imf)
     !thr = MAXVAL(em(1:4)) * 1e-3_dp
     thr = MAX(1e-3_dp,  MAXVAL(em(1:4)) * 1e-2_dp)
     !
     IF(ABS(em(4) - em(3)) < thr ) THEN
        IF(ABS(em(4) - em(2)) < thr ) THEN
           IF(ABS(em(4) - em(1)) < thr ) THEN
              !
              ! em(4) = em(3) = em(2) = em(1)
              !
              w2(1,4) = 0.25_dp * em(4) / ((1.0_dp + em(4)**2))
              w2(2,4) = 0.25_dp         / ((1.0_dp + em(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           ELSE
              !
              ! em(4) = em(3) = em(2)
              !
              w2(1:2,4) = lindhard_1211(em(4),em(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = lindhard_1222(em(1),em(4))
              !
              IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
                 WRITE(*,'(100e15.5)') em(1:4)
                 WRITE(*,'(2e15.5)') w2(1:2,1:4)
                 CALL errore("lindhard", "Stop in lindhard. weighting 4=3=2. imf = ", imf)
              END IF
              !
           END IF
        ELSE IF(ABS(em(2) - em(1)) < thr ) THEN
           !
           ! em(4) = em(3), em(2) = em(1)
           !
           w2(1:2,4) = lindhard_1221(em(4),em(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1221(em(2),em(4))
           w2(1:2,1) = w2(1:2,2)
           !
           IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') em(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              CALL errore("lindhard", "Stop in lindhard. weighting 4=3, 2=1. imf = ", imf)
           END IF
           !
        ELSE
           !
           ! em(4) = em(3)
           !
           w2(1:2,4) = lindhard_1231(em(4),em(1),em(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1233(em(2),em(1),em(4))
           w2(1:2,1) = lindhard_1233(em(1),em(2),em(4))
           !
           IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') em(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              CALL errore("lindhard", "Stop in lindhard. weighting 4=3. imf = ", imf)
           END IF
           !
        END IF
     ELSE IF(ABS(em(3) - em(2)) < thr) THEN
        IF(ABS(em(3) - em(1)) < thr) THEN
           !
           ! em(3) = em(2) = em(1)
           !
           w2(1:2,4) = lindhard_1222(em(4),em(3))
           w2(1:2,3) = lindhard_1211(em(3),em(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') em(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              CALL errore("lindhard", "Stop in lindhard. weighting 3=2=1. imf = ", imf)
           END IF
           !
        ELSE
           !
           ! em(3) = em(2)
           !
           w2(1:2,4) = lindhard_1233(em(4),em(1),em(3))
           w2(1:2,3) = lindhard_1231(em(3),em(1),em(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = lindhard_1233(em(1),em(4),em(3))
           !
           IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') em(1:4)
              WRITE(*,'(2e15.5)') w2(1:2,1:4)
              CALL errore("lindhard", "Stop in lindhard. weighting 3=2. imf = ", imf)
           END IF
           !
        END IF
     ELSE IF(ABS(em(2) - em(1)) < thr) THEN
        !
        ! em(2) = em(1)
        !
        w2(1:2,4) = lindhard_1233(em(4),em(3),em(2))
        w2(1:2,3) = lindhard_1233(em(3),em(4),em(2))
        w2(1:2,2) = lindhard_1231(em(2),em(3),em(4))
        w2(1:2,1) = w2(1:2,2)
        !
        IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') em(1:4)
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           CALL errore("lindhard", "Stop in lindhard. weighting 2=1. imf = ", imf)
        END IF
        !
     ELSE
        !
        ! Different each other.
        !
        w2(1:2,4) = lindhard_1234(em(4),em(1),em(2),em(3))
        w2(1:2,3) = lindhard_1234(em(3),em(1),em(2),em(4))
        w2(1:2,2) = lindhard_1234(em(2),em(1),em(3),em(4))
        w2(1:2,1) = lindhard_1234(em(1),em(2),em(3),em(4))
        !      
        IF(ANY(w2(1:2,1:4) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') em(1:4)
           WRITE(*,'(2e15.5)') w2(1:2,1:4)
           CALL errore("lindhard", "Stop in lindhard. weighting each other. imf = ", imf)
        END IF
        !
     END IF
     !
     w(imf,itetra(1:4)) = CMPLX(w2(1,1:4) /    mf(imf), &
     &                          w2(2,1:4) / (- mf(imf)), dp)
     !
  END DO ! imf
  !
  w(0:nmf,1:4) = - w(0:nmf,1:4)
  !
END SUBROUTINE lindhard
!
! 1, Different each other
!
FUNCTION lindhard_1234(g1,g2,g3,g4) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2, g3, g4
  REAL(dp) :: w(2)
  !
  REAL(dp) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2.0_dp*(3.0_dp*g2**2 - 1.0_dp)*(ATAN(g2) - ATAN(g1)) + (g2**2 - &
  &      3.0_dp)*g2*LOG((1.0_dp + g2**2)/( 1.0_dp + g1**2))
  w2 = -2.0_dp*(g2**2 - 1.0_dp) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2.0_dp*(3.0_dp*g3**2 - 1.0_dp)*(ATAN(g3) - ATAN(g1)) + (g3**2 -  &
  &      3.0_dp)*g3*LOG((1.0_dp + g3**2)/( 1.0_dp + g1**2))
  w3 = -2.0_dp*(g3**2 - 1.0_dp) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2.0_dp*(3.0_dp*g4**2 - 1.0_dp)*(ATAN(g4) - ATAN(g1)) + (g4**2 -  &
  &      3.0_dp)*g4*LOG((1.0_dp + g4**2)/( 1.0_dp + g1**2))
  w4 = -2.0_dp*(g4**2 - 1.0_dp) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2.0_dp*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2.0_dp*(3.0_dp - g2**2)* &
  &    g2*(ATAN(g2) - ATAN(g1)) + (3.0_dp*g2**2 - 1.0_dp)* &
  &    LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w2 = 4.0_dp*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2.0_dp*(3.0_dp - g3**2)* &
  &    g3*(ATAN(g3) - ATAN(g1)) + (3.0_dp*g3**2 - 1.0_dp)* &
  &    LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = 4.0_dp*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2.0_dp*(3.0_dp - g4**2)* &
  &    g4*(ATAN(g4) - ATAN(g1)) + (3.0_dp*g4**2 - 1.0_dp)* &
  &    LOG((1.0_dp + g4**2)/(1.0_dp + g1**2))
  w4 = 4.0_dp*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2.0_dp*(g4 - g2))
  !
END FUNCTION lindhard_1234
!
! 2, g4 = g1
!
FUNCTION lindhard_1231(g1,g2,g3) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2, g3
  REAL(dp) :: w(2)
  !
  REAL(dp) :: w2, w3
  !
  ! Real
  !
  w2 = 2.0_dp*(-1.0_dp + 3.0_dp*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(-3.0_dp + g2**2)*LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w2 = 2.0_dp*(1.0_dp - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2.0_dp*(-1.0_dp + 3.0_dp*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(-3.0_dp + g3**2)*LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = 2.0_dp*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2.0_dp*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2.0_dp* &
  &    g2*(3.0_dp - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1.0_dp + 3.0_dp*g2**2)* &
  &    LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w2 = 4.0_dp*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2.0_dp* &
  &    g3*(3.0_dp - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1.0_dp + 3.0_dp*g3**2)* &
  &    LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = 4.0_dp*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2.0_dp*(g3 - g2))
  !
END FUNCTION lindhard_1231
!
! 3, g4 = g3
!
FUNCTION lindhard_1233(g1, g2, g3) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2, g3
  REAL(dp) :: w(2)
  !
  REAL(dp) :: w2, w3
  !
  ! Real
  !
  w2 = 2.0_dp*(1.0_dp - 3.0_dp*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(3.0_dp - g2**2)*LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w2 = 2.0_dp*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2.0_dp*(1.0_dp - 3.0_dp*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(3.0_dp - g3**2)*LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = 2.0_dp*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4.0_dp*(1.0_dp - 3.0_dp*g1*g3)*(ATAN(g3) - ATAN(g1)) + (3.0_dp*g1 +  &
  &      3.0_dp*g3 - 3.0_dp*g1*g3**2 + g3**3) * LOG((1.0_dp + g3**2)/( &
  &     1.0_dp + g1**2))
  w3 = -4.0_dp*(1.0_dp - g1**2) + w3/(g3 - g1)
  w3 = 4.0_dp*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2.0_dp*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2.0_dp* &
  &    g2*(3.0_dp - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1.0_dp + 3.0_dp*g2**2)* &
  &    LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w2 = 4.0_dp*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2.0_dp* &
  &    g3*(3.0_dp - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1.0_dp + 3.0_dp*g3**2)* &
  &    LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = 4.0_dp*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3.0_dp*g1 - 3.0_dp*g1*g3**2 + 3.0_dp*g3 + g3**3)*(ATAN(g3) -  &
  &      ATAN(g1)) + (3.0_dp*g1*g3 - 1.0_dp)* &
  &    LOG((1.0_dp + g3**2)/(1.0_dp + g1**2))
  w3 = w3/(g3 - g1) - 4.0_dp*g1
  w3 = w3/(g3 - g1) - 2.0_dp
  w3 = (2.0_dp*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2.0_dp*(g3 - g2))
  !
END FUNCTION lindhard_1233
!
! 4, g4 = g1 and g3 = g2
!
FUNCTION lindhard_1221(g1,g2) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2
  REAL(dp) :: w(2)
  !
  ! Real
  !
  w(1) = -2.0_dp*(-1.0_dp + 2.0_dp*g1*g2 + g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (g1 + 2.0_dp*g2 - g1*g2**2)* &
  &    LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w(1) = 2.0_dp*(-1.0_dp + g1**2) + w(1)/(g2 - g1)
  w(1) = 3.0_dp*g1 + w(1)/(g2 - g1)
  w(1) = 2.0_dp + (3.0_dp*w(1))/(g2 - g1)
  w(1) = w(1)/(2.0_dp*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2.0_dp*(g1 + 2.0_dp*g2 - g1*g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-1.0_dp + 2.0_dp*g1*g2 + g2**2)* &
  &    LOG((1 + g2**2)/(1 + g1**2))
  w(2) = -4.0_dp*g1 + w(2)/(g2 - g1)
  w(2) = -3.0_dp + w(2)/(g2 - g1)
  w(2) = (3.0_dp*w(2))/(2.0_dp*(g2 - g1)**2)
  !
END FUNCTION lindhard_1221
!
! 5, g4 = g3 = g2
!
FUNCTION lindhard_1222(g1,g2) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2
  REAL(dp) :: w(2)
  !
  ! Real
  !
  w(1) = 2.0_dp*(-1.0_dp + g1**2 + 2.0_dp*g1*g2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-2.0_dp*g1 - g2 + g1**2*g2) * LOG((1.0_dp + g2**2)/( &
  &     1.0_dp + g1**2))
  w(1) = 2.0_dp*(1.0_dp - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1.0_dp - (3.0_dp*w(1))/(g2 - g1)
  w(1) = w(1)/(2.0_dp*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2.0_dp*(-2.0_dp*g1 - g2 + g1**2*g2)*(ATAN(g2) - ATAN(g1)) + (1.0_dp - &
  &       g1**2 - 2.0_dp*g1*g2) * LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w(2) = 4.0_dp*g1 + w(2)/(g2 - g1)
  w(2) = 1.0_dp + w(2)/(g2 - g1)
  w(2) = (3.0_dp*w(2))/(2.0_dp*(g2 - g1)**2)
  !
END FUNCTION lindhard_1222
!
! 6, g4 = g3 = g1
!
FUNCTION lindhard_1211(g1,g2) RESULT(w)
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in) :: g1, g2
  REAL(dp) :: w(2)
  !
  ! Real
  !
  w(1) = 2.0_dp*(3.0_dp*g2**2 - 1.0_dp)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(g2**2 - 3.0_dp)*LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w(1) = 2.0_dp*(1.0_dp - g1**2) + w(1)/(g2 - g1)
  w(1) = -5.0_dp*g1 + w(1)/(g2 - g1)
  w(1) = -11.0_dp + (3.0_dp*w(1))/(g2 - g1)
  w(1) = w(1)/(6.0_dp*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2.0_dp*g2*(-3.0_dp + g2**2)*(ATAN(g2) - ATAN(g1)) + (1.0_dp -  &
  &      3.0_dp*g2**2)*LOG((1.0_dp + g2**2)/(1.0_dp + g1**2))
  w(2) = 4.0_dp*g2 + w(2)/(g2 - g1)
  w(2) = 1.0_dp + w(2)/(g2 - g1)
  w(2) = w(2)/(2.0_dp*(g2 - g1)**2)
  !
END FUNCTION lindhard_1211
!
END MODULE sctk_tetra
