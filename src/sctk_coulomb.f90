!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_coulomb
  !
  IMPLICIT NONE
  !
CONTAINS
!>
!> Allocate Kel & ikq
!>
SUBROUTINE alloc_Kel()
  !
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE sctk_val, ONLY : gindx, gq2, Kel, nftot, nqbz, nmf, nci, lsf, mf
  !
  IMPLICIT NONE
  !
  INTEGER :: imf
  REAL(dp) :: x0
  !
  nmf = 2*nci - 3
  !
  ALLOCATE(gq2(nftot), gindx(nftot), mf(nmf))
  IF(lsf>0) THEN
     WRITE(stdout,'(7x,"RAM for Kel per process : ",e10.2," GB")') &
     &     REAL((nmf+2)*2,dp)*REAL(nbnd,dp)**2*REAL(nqbz,dp)*8.0e-9_dp
     ALLOCATE(Kel(0:nmf+1,nbnd,nbnd,nqbz,2))
  ELSE
     WRITE(stdout,'(7x,"RAM for Kel per process : ",e10.2," GB")') &
     &     REAL((nmf+2),dp)*REAL(nbnd,dp)**2*REAL(nqbz,dp)*8.0e-9_dp
     ALLOCATE(Kel(0:nmf+1,nbnd,nbnd,nqbz,1))
  END IF
  !
  x0 = COS(pi / REAL(2 * nci, dp))
  !
  WRITE(stdout,'(7x,"Index, Frequency [Ry]")')
  WRITE(stdout,'(9x,i0,5x,e15.7)') 0, 0.0_dp
  DO imf = 1, nmf
     mf(imf) = - COS(REAL(imf+1, dp) * pi / REAL(2 * nci, dp))
     mf(imf) = (x0 + mf(imf)) / (x0 - mf(imf))
     WRITE(stdout,'(9x,i0,5x,e15.7)') imf, mf(imf)
  END DO
  WRITE(stdout,'(9x,i0,5x,"Infinity")') nmf+1  
  !
END SUBROUTINE alloc_Kel
!>
!>
!>
SUBROUTINE circular_shift_wrapper(wfc)
  !
  USE kinds, ONLY : DP
  USE sctk_val, ONLY : nkpe, nftot
  USE wvfct, ONLY : nbnd
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_circular_shift_left
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(IN) :: wfc(nftot*npol*nbnd,nkpe)
  !
  CALL mp_circular_shift_left( wfc, 1, world_comm )
  !
END SUBROUTINE circular_shift_wrapper
!>
!> Screened coulomb interaction
!>
SUBROUTINE prepare_q()
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE cell_base, ONLY : bg, tpiba
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout
  USE gvecw, ONLY : ecutwfc
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE cell_base, ONLY : at
  USE control_ph, ONLY : current_iq
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : gindx, gq2, nf, nftot, ngv, nqbz, nkpe, &
  &                    wfc1, wfc1q, wfc2, wfc2q, wscr, &
  &                    ngv0, ngv1, lsf, nmf
  ! 
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp_full
  IMPLICIT NONE
  !
  INTEGER :: ik, jk, ib, i1, i2, i3, ikv(3), jkv(3), g0(3), ir(3), ifft, ig, igv(3), &
  &          cnt(0:nproc - 1), dsp(0:nproc - 1), jkindx(nqbz), org, ipe, iqv(3), ipol
  REAL(dp) :: gv(3), theta, gq20
  COMPLEX(dp) :: phase(nftot), wfctmp(nftot,npol,nbnd,nkpe)
  !
  ! |G+q|^2
  !
  ifft = 0
  ngv = 0
  DO i3 = 1, nf(3)
     DO i2 = 1, nf(2)
        DO i1 = 1, nf(1)
           !
           ifft = ifft + 1
           !
           igv(1:3) = (/i1, i2, i3/) - 1
           WHERE(igv(1:3)*2 >= nf(1:3)) igv(1:3) = igv(1:3) - nf(1:3)
           gv(1:3) = MATMUL(bg(1:3,1:3), REAL(igv(1:3), dp)) * tpiba
           gv(1:3) = gv(1:3) - x_q(1:3, current_iq) * tpiba
           !
           gq20 = DOT_PRODUCT(gv(1:3), gv(1:3))
           !
           IF(ecutwfc < 1e-10_dp .OR. gq20 < ecutwfc) THEN
              !
              ngv = ngv + 1
              gq2(ngv) = gq20 / (8.0_dp * pi)
              gindx(ngv) = ifft
              !
           END IF
           !
        END DO ! i1 = 1, nf(1)
     END DO ! i2 = 1, nf(2)
  END DO ! i3 = 1, nf(3)
  !
  WRITE(stdout,'(9x,"Number of PWs for W : ",i0)') ngv
  !
  CALL divide(world_comm, ngv, ngv0, ngv1)
  IF(lsf >0) THEN
     WRITE(stdout,'(9x,"Total RAM for Vscr : ",e10.2," GB")') &
     &     REAL(ngv,dp)**2*REAL((nmf+1)*2*npol,dp)*16.0e-9_dp
     ALLOCATE(wscr(ngv, ngv0:ngv1, 0:nmf, 2*npol))
  ELSE
     WRITE(stdout,'(9x,"Total RAM for Vscr : ",e10.2," GB")') &
     &     REAL(ngv,dp)**2*REAL((nmf+1),dp)*16.0e-9_dp
     ALLOCATE(wscr(ngv, ngv0:ngv1, 0:nmf, 1))
  END IF
  !
  ! Prepare wave functions with pahse shift
  !
  CALL cnt_and_dsp_full(nqbz, cnt, dsp)
  !
  iqv(1:3) = NINT(MATMUL(x_q(1:3,current_iq), at(1:3,1:3)) * REAL((/nq1, nq2, nq3/), dp) - 0.5_dp)
  !
  DO ik = 1, nqbz
     !
     ikv(1) = (ik - 1) / (nq3*nq2)
     ikv(2) = (ik - 1 - ikv(1)*nq2*nq3) / nq3
     ikv(3) =  ik - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
     !
     WHERE(ikv(1:3)*2 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
     !
     jkv(1:3) = ikv(1:3) + iqv(1:3)
     jkv(1:3) = MODULO(jkv(1:3), (/nq1,nq2,nq3/))
     WHERE(jkv(1:3)*2 + 1 >= (/nq1,nq2,nq3/)) jkv(1:3) = jkv(1:3) - (/nq1,nq2,nq3/)
     !
     g0(1:3) = (jkv(1:3) - ikv(1:3) - iqv(1:3)) / (/nq1,nq2,nq3/)
     !
     jkv(1:3) = MODULO(jkv(1:3), (/nq1,nq2,nq3/))
     jk = 1 + jkv(3) + jkv(2)*nq3 + jkv(1)*nq2*nq3
     jkindx(ik) = jk
     !
     IF(ik <= dsp(mpime) .OR. dsp(mpime) + cnt(mpime) < ik) CYCLE
     !
     ig = 0
     DO i3 = 1, nf(3)
        DO i2 = 1, nf(2)
           DO i1 = 1, nf(1)
              !
              ig = ig + 1
              !
              ir(1:3) = (/i1, i2, i3/) - 1
              ir(1:3) = ir(1:3) * g0(1:3)
              !             
              theta = SUM(REAL(ir(1:3), dp) / REAL(nf(1:3), dp))
              theta = - 2.0_dp * pi * theta
              !
              phase(ig) = CMPLX(COS(theta), SIN(theta), KIND=dp)
              !
           END DO ! i1
        END DO ! i2
     END DO ! i3
     !
     DO ib = 1, nbnd
        DO ipol = 1, npol
           !
           wfc1q(       1:nftot,ipol,ib,ik - dsp(mpime)) = &
           & wfc1(      1:nftot,ipol,ib,ik - dsp(mpime)) * phase(1:nftot)
           wfc2q(       1:nftot,ipol,ib,ik - dsp(mpime)) = &
           & CONJG(wfc2(1:nftot,ipol,ib,ik - dsp(mpime)))
           !
        END DO ! ipol
     END DO ! ib
     !
  END DO ! ik
  !
  !
  !
  wfctmp(1:nftot,1:npol,1:nbnd,1:nkpe) = wfc2q(1:nftot,1:npol,1:nbnd,1:nkpe)
  wfc2q( 1:nftot,1:npol,1:nbnd,1:nkpe) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  DO ipe = 1, nproc
     !
     CALL circular_shift_wrapper(wfctmp)
     !
     org = MODULO(mpime + ipe, nproc)
     !
     DO ik = 1, cnt(mpime)
        !
        IF(jkindx(dsp(mpime) + ik) <= dsp(org) .OR. &
        &  dsp(org) + cnt(org) < jkindx(dsp(mpime) + ik)) CYCLE
        !
        wfc2q(   1:nftot,1:npol,1:nbnd,ik) = &
        & wfctmp(1:nftot,1:npol,1:nbnd,jkindx(dsp(mpime) + ik) - dsp(org))
        !
     END DO
     !
  END DO
  !
END SUBROUTINE prepare_q
!>
!> Calculation of weight function 
!> @f$\frac{f(1-f')}{\varepsilon - \evarepsilon'}@f$
!>
SUBROUTINE fermi_factor(wght)
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE ktetra, ONLY : ntetra, tetra
  USE cell_base, ONLY : at
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE start_k, ONLY : nk1, nk2, nk3
  USE sctk_val, ONLY : nqbz, nmf
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE wvfct, ONLY : et, nbnd
  USE klist, ONLY : nks
  USE control_ph, ONLY : current_iq
  !
  USE sctk_tetra, ONLY : tetraweight, interpol_indx
  !
  IMPLICIT NONE
  !
  COMPLEX(dp),INTENT(OUT) :: wght((nmf+1)*nbnd*nbnd,nqbz) !< integration weight
  !
  INTEGER :: ntetra0, ntetra1, nks0, it, ik, ii, dik(3), ikv(3), &
  &          indx2(20, ntetra), indx3(20 * ntetra), kintp(8), ikq
  REAL(dp) :: kv(3), wintp(1,8)
  COMPLEX(dp),ALLOCATABLE :: wghtd(:,:,:)
  !
  dik(1:3) = NINT(MATMUL(x_q(1:3,current_iq), at(1:3,1:3)) * REAL((/nk1, nk2, nk3/), dp))
  !
  DO ik = 1, nks
     !
     ikv(1) = (ik - 1) / (nk3*nk2)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     !
     ikv(1:3) = ikv(1:3) + dik(1:3)
     ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/))
     ikq = 1 + ikv(3) + ikv(2)*nk3 + ikv(1)*nk2*nk3
     !
     et(1:nbnd, nks+ik) = et(1:nbnd, ikq)
     !
  END DO
  !
  indx2(1:20,1:ntetra) = 0
  indx3(1:20 * ntetra) = 0
  !
  CALL divide(world_comm, ntetra,ntetra0,ntetra1)
  !
  nks0 = 0
  DO it = ntetra0, ntetra1
     !
     DO ii = 1, 20
        !
        DO ik = 1, nks0
           !
           IF(tetra(ii,it) == indx3(ik)) THEN
              !
              indx2(ii,it) = ik
              GOTO 10
              !
           END IF
           !
        END DO
        !
        nks0 = nks0 + 1
        indx2(ii,it) = nks0
        indx3(nks0) = tetra(ii,it)
        !
10      CONTINUE
        !
     END DO
     !
  END DO
  !
  ALLOCATE(wghtd((nmf+1)*nbnd*nbnd,1,nks0))
  !
  CALL tetraweight(nks0,indx2,wghtd)
  !
  ! Interpolation of weight
  !
  wght(1:(nmf+1)*nbnd*nbnd,1:nqbz) = 0.0_dp
  DO ik = 1, nks0
     !
     ikv(1) = (indx3(ik) - 1) / (nk3*nk2)
     ikv(2) = (indx3(ik) - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  indx3(ik) - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     !
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
     wght(1:(nmf+1)*nbnd*nbnd,kintp(1:8)) = wght(1:(nmf+1)*nbnd*nbnd,             kintp(1:8)) &
     &                            + MATMUL(wghtd(1:(nmf+1)*nbnd*nbnd,1:1,ik), wintp(1:1,1:8))
  END DO
  !
  CALL mp_sum( wght, world_comm )
  !
  DEALLOCATE(wghtd)
  !
END SUBROUTINE fermi_factor
!>
!> Compute screened interaction
!>
SUBROUTINE make_scrn()
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime, nproc
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : gindx, nf, nftot, ngv, nqbz, nmf, &
  &                    wfc1q, wfc2q, wscr, lsf, ngv0, ngv1
  !
  USE sctk_cnt_dsp, ONLY :  cnt_and_dsp_full
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jb, imf, org, ipe, ipol, npol2, nwscr, &
  &          cnt(0:nproc - 1), dsp(0:nproc - 1)
  !
  COMPLEX(dp) :: wght(0:nmf,nbnd,nbnd,nqbz), rho0(nftot, 4)
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:,:)
  !
  CALL start_clock("make_scrn")
  !
  CALL cnt_and_dsp_full(nqbz, cnt, dsp)
  !
  IF(lsf>0) THEN
     IF(npol == 2) THEN
        npol2 = 4
        nwscr = 4
     ELSE
        npol2 = 1
        nwscr = 2
     END IF
  ELSE
     npol2 = 1
     nwscr = 1
  END IF
  !
  ALLOCATE(rho1(ngv,nbnd,npol2), rho2(ngv0:ngv1,nbnd,npol2))
  wscr(1:ngv, ngv0:ngv1, 0:nmf, 1:nwscr) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Calc f * (1 - f') / (e - e' + iw)
  !
  CALL fermi_factor(wght)
  !
  ! Calc. Chi
  !
  DO ipe = 1, nproc
     !
     CALL circular_shift_wrapper(wfc1q)
     CALL circular_shift_wrapper(wfc2q)
     !
     org = MODULO(mpime + ipe, nproc)
     !
     DO ik = 1, cnt(org)
        !
        DO ib = 1, nbnd
           !
           DO jb = 1, nbnd
              !
              rho0(1:nftot, 1:npol2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
              DO ipol = 1, npol
                 !
                 rho0(1:nftot,1) = rho0( 1:nftot,1) &
                 &               + wfc1q(1:nftot,ipol,ib,ik) &
                 &               * wfc2q(1:nftot,ipol,jb,ik)
                 !
                 IF(npol2 == 4) THEN
                    !
                    rho0(1:nftot,2) = rho0( 1:nftot,2) &
                    &               + wfc1q(1:nftot,  ipol,ib,ik) &
                    &               * wfc2q(1:nftot,3-ipol,jb,ik)
                    !
                    rho0(1:nftot,3) = rho0( 1:nftot,3) &
                    &               + REAL(3-2*ipol, dp) &
                    &               * wfc1q(1:nftot,  ipol,ib,ik) &
                    &               * wfc2q(1:nftot,3-ipol,jb,ik)
                    !
                    rho0(1:nftot,4) = rho0( 1:nftot,4) &
                    &               + REAL(3-2*ipol, dp) &
                    &               * wfc1q(1:nftot,ipol,ib,ik) &
                    &               * wfc2q(1:nftot,ipol,jb,ik)
                    !
                 END IF
                 !
              END DO
              !
              DO ipol = 1, npol2
                 CALL cfft3d(rho0(1:nftot, ipol), &
                 & nf(1), nf(2), nf(3), nf(1), nf(2), nf(3), 1, -1)
                 !
                 rho1(1:ngv,jb,ipol) = CONJG(rho0(gindx(1:ngv), ipol))
              END DO
              !
           END DO ! jb
           !
           DO imf = 0, nmf
              !
              DO jb = 1, nbnd
                 rho2(ngv0:ngv1,jb,1:npol2) = REAL(wght(imf,jb,ib,dsp(org) + ik),dp) &
                 &                          * CONJG(rho1(ngv0:ngv1,jb,1:npol2))
              END DO ! jb = 1, nbnd
              !
              DO ipol = 1, npol2
                 CALL zgemm("N", "T", ngv, ngv1-ngv0+1, nbnd, &
                 &  (1.0_dp, 0.0_dp), rho1(1:ngv,            1:nbnd, ipol), ngv, &
                 &                    rho2(       ngv0:ngv1, 1:nbnd, ipol), ngv1-ngv0+1, &
                 &  (1.0_dp, 0.0_dp), wscr(1:ngv, ngv0:ngv1,    imf, ipol), ngv  )
              END DO
              !
           END DO ! imf = 0, nmf
           !
        END DO ! ib = 1, nbnd
        !
     END DO ! ik = 1, cnt(org)
     !
  END DO ! ipe = 1, nproc
  !
  IF(nwscr == 2) &
  &  wscr(1:ngv, ngv0:ngv1, 0:nmf, 2) = wscr(1:ngv, ngv0:ngv1, 0:nmf, 1)
  DEALLOCATE(rho1, rho2)
  !
  CALL stop_clock("make_scrn")
  !
END SUBROUTINE make_scrn
!>
!> Calc. Kel
!>
SUBROUTINE make_Kel()
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE mp_world, ONLY : mpime, nproc
  USE mp, ONLY : mp_sum
  USE fft_scalar, ONLY : cfft3d
  USE mp_world, ONLY : world_comm
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : gindx, gq2, nf, nftot, ngv, nqbz, nmf, &
  &                    wfc1q, wfc2q, ngv0, ngv1, Kel, Wscr, lsf
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp_full
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jb, imf, org, ipe, ipol, jpol, npol2, nwscr, &
  &          cnt(0:nproc - 1), dsp(0:nproc - 1)
  COMPLEX(dp) :: rho0(nftot,4), Kel0(4)
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:,:), rho3(:,:,:)
  COMPLEX(DP),EXTERNAL :: zdotc
  !
  CALL start_clock("make_Kel")
  !
  CALL cnt_and_dsp_full(nqbz, cnt, dsp)
  !
  IF(lsf>0) THEN
     IF(npol == 2) THEN
        npol2 = 4
        nwscr = 4
     ELSE
        npol2 = 1
        nwscr = 2
     END IF
     Kel(0:nmf + 1,1:nbnd,1:nbnd,1:nqbz,1:2) = 0.0_dp
  ELSE
     npol2 = 1
     nwscr = 1
     Kel(0:nmf + 1,1:nbnd,1:nbnd,1:nqbz,1) = 0.0_dp
  END IF
  !
  ALLOCATE(rho1(ngv,nbnd,npol2), rho2(ngv0:ngv1,nbnd,npol2), rho3(ngv,nbnd,nwscr))
  !
  DO ipe = 1, nproc
     !
     CALL circular_shift_wrapper(wfc1q)
     CALL circular_shift_wrapper(wfc2q)
     !
     org = MODULO(mpime + ipe, nproc)
     !
     DO ik = 1, cnt(org)
        !
        DO ib = 1, nbnd
           !
           DO jb = 1, nbnd
              !
              rho0(1:nftot, 1:npol2) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
              DO ipol = 1, npol
                 !
                 rho0(1:nftot,1) = rho0( 1:nftot,1) &
                 &               + wfc1q(1:nftot,ipol,ib,ik) &
                 &               * wfc2q(1:nftot,ipol,jb,ik)
                 !
                 IF(npol2 == 4) THEN
                    !
                    rho0(1:nftot,2) = rho0( 1:nftot,2) &
                    &               + wfc1q(1:nftot,  ipol,ib,ik) &
                    &               * wfc2q(1:nftot,3-ipol,jb,ik)
                    !
                    rho0(1:nftot,3) = rho0( 1:nftot,3) &
                    &               + REAL(3-2*ipol, dp) &
                    &               * wfc1q(1:nftot,  ipol,ib,ik) &
                    &               * wfc2q(1:nftot,3-ipol,jb,ik)
                    !
                    rho0(1:nftot,4) = rho0( 1:nftot,4) &
                    &               + REAL(3-2*ipol, dp) &
                    &               * wfc1q(1:nftot,ipol,ib,ik) &
                    &               * wfc2q(1:nftot,ipol,jb,ik)
                    !
                 END IF ! (npol2 == 4)
                 !
              END DO
              !
              DO ipol = 1, npol2
                 CALL cfft3d (rho0(1:nftot,ipol), &
                 & nf(1), nf(2), nf(3), nf(1), nf(2), nf(3), 1, -1)
                 !
                 rho1(   1:ngv, jb,ipol) = rho0(gindx(1:ngv),ipol)
                 rho2(ngv0:ngv1,jb,ipol) = rho1(ngv0:ngv1,jb,ipol)
              END DO
              !
           END DO ! jb = 1, nbnd
           !
           DO imf = 0, nmf
              !
              DO ipol = 1, nwscr
                 IF(nwscr == 2) THEN
                    jpol = 1
                 ELSE
                    jpol = ipol
                 END IF
                 CALL zgemm("N", "N", ngv, nbnd, ngv1-ngv0+1, &
                 &          (1.0_dp, 0.0_dp),Wscr(1:ngv,ngv0:ngv1,   imf,ipol),ngv, &
                 &                           rho2(      ngv0:ngv1,1:nbnd,jpol),ngv1-ngv0+1, &
                 &          (0.0_dp, 0.0_dp),rho3(1:ngv,          1:nbnd,ipol),ngv  )
              END DO ! ipol = 1, nwscr
              !
              DO jb = 1, nbnd
                 !
                 DO ipol = 1, nwscr
                    IF(nwscr == 2) THEN
                       jpol = 1
                    ELSE
                       jpol = ipol
                    END IF
                    Kel0(ipol) = zdotc(ngv, rho1(1:ngv,jb,jpol), 1, &
                    &                       rho3(1:ngv,jb,ipol), 1)
                 END DO
                 !
                 Kel(imf,jb,ib,dsp(org) + ik,1) = REAL(    Kel0(1),    dp) / omega
                 IF(nwscr == 2) THEN
                    Kel(imf,jb,ib,dsp(org) + ik,2) = REAL(    Kel0(2),    dp) / omega
                 ELSE IF(nwscr == 4) THEN
                    Kel(imf,jb,ib,dsp(org) + ik,2) = REAL(SUM(Kel0(2:4)), dp) / omega
                 END IF
                 !
              END DO ! jb = 1, nbnd
              !
           END DO ! imf = 0, nmf
           !
           ! Infinite frequency -> Bare Coulomb
           !
           DO jb = 1, nbnd
              Kel(nmf + 1,jb,ib,dsp(org) + ik,1) = SUM(REAL(rho2(ngv0:ngv1,jb,1) &
              &                                     * CONJG(rho2(ngv0:ngv1,jb,1)), dp) &
              &                                     /        gq2(ngv0:ngv1))  / omega
           END DO ! jb = 1, nbnd
           !
        END DO ! ib = 1, nbnd
        !
     END DO ! ik = 1, cnt(org)
     !
  END DO ! ipe = 1, nproc
  !
  CALL mp_sum( Kel, world_comm )
  !
  DEALLOCATE(rho1, rho2, rho3)
  !
  CALL stop_clock("make_Kel")
  !
END SUBROUTINE make_Kel
!>
!>
!>
SUBROUTINE Coulomb_parameter()
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : nbnd, et
  USE lsda_mod,   ONLY : nspin
  USE ener, ONLY : ef
  USE ktetra,     ONLY : opt_tetra_dos_t
  USE klist,      ONLY : nks, nkstot
  USE io_global, ONLY : stdout
  USE el_phon, ONLY : elph_nbnd_min
  USE mp_world, ONLY : world_comm
  USE fermisurfer_common, ONLY : b_low, b_high
  USE ktetra, ONLY : ntetra
  USE sctk_tetra, ONLY : interpol_indx
  USE start_k, ONLY : nk1, nk2, nk3
  USE disp,  ONLY : nq1, nq2, nq3
  USE elph_tetra_mod, ONLY : elph_tetra_delta1
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : nmf, mf, Kel, lsf
  !
  IMPLICIT NONE
  !
  INTEGER :: ntetra0, ntetra1, ii, ikv(3), ik, imf, kintp(8), nkel, ikel
  REAL(dp) :: dost(2), kv(3), mu(0:nmf+1,2), wintp(8), scale_mu
  REAL(dp),allocatable :: wght(:,:,:)
  !
  WRITE(stdout,'(/,9x,"Averaged potential (mu)",/)')
  !
  IF(lsf>0) THEN
     nkel = 2
  ELSE
     nkel = 1
  END IF
  !
  ALLOCATE(wght(b_low:b_high, b_low:b_high, nkstot))
  !
  ! DOS
  !
  CALL opt_tetra_dos_t (et, nspin, nbnd, nks, ef, dost)
  WRITE(stdout,'(11x,"DOS(EF) = ",e14.4," [/Ry]")') dost(1)
  !
  ! Double delta
  !
  CALL divide(world_comm, ntetra,ntetra0,ntetra1)
  !
  elph_nbnd_min = b_low
  CALL elph_tetra_delta1(b_high-b_low+1,nks,ntetra0,ntetra1,et,wght)
  !
  ! Interpolation
  !
  mu(0:nmf+1, 1:nkel) = 0.0_dp
  DO ik = 1, nks
     !
     ikv(1) = (ik - 1) / (nk3*nk2)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     !
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
     !
     DO ikel = 1, nkel
        DO imf = 0, nmf + 1
           DO ii = 1, 8
              mu(imf,ikel) = mu(imf,ikel) + &
              & SUM(Kel(imf,b_low:b_high,b_low:b_high,       kintp(ii), ikel) &
              &     * wght( b_low:b_high,b_low:b_high,ik)) * wintp(ii)
           END DO ! ii = 1, 8
        END DO ! imf = 0, nmf + 1
     END DO ! ikel = 1, nkel
     !
  END DO ! ik = 1, nks
  !
  scale_mu = dost(1) / SUM(wght( b_low:b_high,b_low:b_high,1:nks))
  mu(0:nmf+1, 1:nkel) = mu(0:nmf+1, 1:nkel) * scale_mu
  IF(npol == 1) mu(0:nmf+1, 1:nkel) = mu(0:nmf+1, 1:nkel) * 0.5_dp
  !
  DO ikel = 1, nkel
     !
     IF(ikel == 1) THEN
        WRITE(stdout,'(11x,"Frequency [Ry], Coulomb parameter(mu)")')
     ELSE
        WRITE(stdout,'(11x,"Frequency [Ry], Spin-F. parameter(mu)")')
     END IF
     !
     WRITE(stdout,'(13x,"           0",e12.4)') mu(0, ikel)
     DO imf = 1, nmf
        WRITE(stdout,'(13x,e12.4,e12.4)') mf(imf), mu(imf, ikel)
     END DO
     WRITE(stdout,'(13x,"    Infinity",e12.4)') mu(nmf+1, ikel)
     !
  END DO ! ikel = 1, nkel
  !
  DEALLOCATE(wght)
  !
END SUBROUTINE Coulomb_parameter
!>
!> Compute parameter for Chebyshev interpolation
!>
SUBROUTINE chebyshev_interpol()
  !
  USE constants, ONLY : pi
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum, mp_max
  USE io_global, ONLY : stdout
  !
  USE sctk_val, ONLY : Kel, nqbz, nmf, lsf, mf, nci
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jb, imf, jmf, nqbz0, nqbz1, ikel, nkel
  REAL(dp) :: coef(nci,nci), Kel0(nci), Kel1(0:nmf+1), cheb(1:nci,0:nmf+1), &
  &           err(0:nmf+1), aveer(0:nmf+1), maxer(0:nmf+1), mf0, x0
  !
  CALL divide(world_comm, nqbz, nqbz0, nqbz1)
  !
  IF(lsf>0) THEN
     nkel = 2
  ELSE
     nkel = 1
  END IF
  !
  Kel(0:nmf + 1,1:nbnd,1:nbnd,      1:nqbz0-1, 1:nkel) = 0.0_dp
  Kel(0:nmf + 1,1:nbnd,1:nbnd,nqbz1+1:nqbz,    1:nkel) = 0.0_dp
  !
  do imf = 1, nci
     do jmf = 1, nci
        coef(jmf,imf) = 2.0_dp / real(nci, dp) &
        &  * COS(REAL((2*(nci-imf) + 1) * (jmf-1), dp) * pi / REAL(2*nci, dp))
     end do ! jmf
  end do ! imf
  coef(1,1:nci) = coef(1,1:nci) * 0.5_dp
  !
  x0 = COS(pi / REAL(2 * nci, dp))
  !
  mf0 = ACOS(-x0)
  DO jmf = 1, nci
     cheb(jmf,0) = COS(REAL(jmf-1,dp) * mf0)
  END DO
  DO imf = 1, nmf
     mf0 = x0 * (mf(imf) - 1.0_dp) / (mf(imf) + 1.0_dp)
     mf0 = ACOS(mf0)
     DO jmf = 1, nci
        cheb(jmf,imf) = COS(REAL(jmf-1,dp) * mf0)
     END DO
  END DO
  mf0 = ACOS(x0)
  DO jmf = 1, nci
     cheb(jmf,nmf+1) = COS(REAL(jmf-1,dp) * mf0)
  END DO
  !
  DO ikel = 1, nkel
     !
     aveer(0:nmf+1) = 0.0_dp
     maxer(0:nmf+1) = 0.0_dp
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(nqbz0,nqbz1,nbnd,nmf,coef,Kel,nci,ikel,cheb,aveer,maxer) &
     !$OMP & PRIVATE(ik,ib,jb,Kel0,Kel1,err)
     !
     !$OMP DO REDUCTION(+:aveer) REDUCTION(MAX:maxer)
     DO ik = nqbz0, nqbz1
        DO ib = 1, nbnd
           DO jb = 1, nbnd
              !
              Kel0(1:nci) = 0.0_dp
              DO imf = 1, nci
                 DO jmf = 1, nci
                    Kel0(imf) = Kel0(imf) + coef(imf,jmf)*Kel((jmf-1)*2,jb,ib,ik,ikel)
                 END DO
              END DO
              !
              Kel1(0:nmf+1) = MATMUL(Kel0(1:nci), cheb(1:nci, 0:nmf+1))
              !
              err(0:nmf+1) = ABS((Kel1(0:nmf+1)  -     Kel(0:nmf+1,jb,ib,ik,ikel))) &
              &  * 2.0_dp / (ABS( Kel1(0:nmf+1)) + ABS(Kel(0:nmf+1,jb,ib,ik,ikel)))
              !
              aveer(0:nmf+1) = aveer(0:nmf+1) + err(0:nmf+1)
              maxer(0:nmf+1) = MAX(maxer(0:nmf+1), err(0:nmf+1))
              !
              Kel(0:nci-1,jb,ib,ik,ikel) = Kel0(1:nci)
              !
           END DO  ! jb
        END DO ! ib
     END DO ! ik
     !$OMP END DO
     !$OMP END PARALLEL
     !
     IF(ikel == 1) THEN
        WRITE(stdout,'(/,9x,"Verify the Chebyshev interpolation for Coulom",/)')
     ELSE
        WRITE(stdout,'(/,9x,"Verify the Chebyshev interpolation for S-F",/)')
     END IF
     WRITE(stdout,'(11x," Frequency   Max. Err.[%]   Ave. Err.[%]")')
     !
     CALL mp_sum( aveer, world_comm )
     CALL mp_max( maxer, world_comm )
     !
     aveer(0:nmf+1) = aveer(0:nmf+1) / REAL(nqbz*nbnd*nbnd, dp)
     maxer(0:nmf+1) = maxer(0:nmf+1)
     !
     WRITE(stdout,'(13x,e15.5,f15.3,f15.3)') 0.0_dp, &
     &     100.0_dp * maxer(0), 100.0_dp * aveer(0)
     DO imf = 1, nmf
        WRITE(stdout,'(13x,e15.5,f15.3,f15.3)') mf(imf), &
        &   100.0_dp * maxer(imf), 100.0_dp * aveer(imf)
     END DO
     IF(ikel == 1) &
     & WRITE(stdout,'(13x,"  Infinity",f15.3,f15.3)') &
     &    100.0_dp * maxer(nmf+1), 100.0_dp * aveer(nmf+1)
     !
  END DO ! ikel
  !
  CALL mp_sum( Kel, world_comm )
  !
END SUBROUTINE chebyshev_interpol
!>
!> Output to file
!>
SUBROUTINE write_Kel()
  !
  USE wvfct, ONLY : nbnd
  USE mp_world, ONLY : mpime
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE cell_base, ONLY : at
  USE control_ph, ONLY : current_iq
  !
  USE sctk_val, ONLY : Kel, nqbz, lsf, nci
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  IF(mpime == 0) THEN
     !
     OPEN(fo, file = "vel"//TRIM(int_to_char(current_iq))//".dat", form = "unformatted")
     !
     WRITE(fo) nq1, nq2, nq3
     WRITE(fo) nbnd
     WRITE(fo) MATMUL(x_q(1:3, current_iq), at(1:3, 1:3))
     WRITE(fo) nci
     !
     WRITE(fo) Kel(0:nci-1,1:nbnd,1:nbnd,1:nqbz,1)
     IF(lsf>0) THEN
        WRITE(fo) Kel(0:nci-1,1:nbnd,1:nbnd,1:nqbz,2)
     END IF
     !
     CLOSE(fo)
     !
  END IF
  !
END SUBROUTINE write_Kel
!
END MODULE sctk_coulomb
