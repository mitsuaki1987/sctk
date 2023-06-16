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
SUBROUTINE Kel_frequency()
  !
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout
  USE kinds, ONLY : dp
  USE sctk_val, ONLY : gindx, gq2, nmf, nci, mf
  USE exx, ONLY : dfftt
  !
  IMPLICIT NONE
  !
  INTEGER :: imf
  REAL(dp) :: x0
  !
  nmf = 2*nci - 3
  !
  ALLOCATE(gq2(dfftt%nnr), gindx(dfftt%nnr), mf(nmf))
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
END SUBROUTINE Kel_frequency
!>
!>
!>
SUBROUTINE circular_shift_wrapper_c(n1, n2, comm, wfc)
  !
  USE kinds, ONLY : DP
  USE mp, ONLY : mp_circular_shift_left
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n1, n2, comm
  COMPLEX(DP),INTENT(INOUT) :: wfc(n1,n2)
  !
  CALL mp_circular_shift_left( wfc, 1, comm )
  !
END SUBROUTINE circular_shift_wrapper_c
!>
!>
!>
SUBROUTINE circular_shift_wrapper_r(n1, n2, comm, wfc)
  !
  USE kinds, ONLY : DP
  USE mp, ONLY : mp_circular_shift_left
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n1, n2, comm
  REAL(DP),INTENT(INOUT) :: wfc(n1,n2)
  !
  CALL mp_circular_shift_left( wfc, 1, comm )
  !
END SUBROUTINE circular_shift_wrapper_r
!>
!> Screened coulomb interaction
!>
SUBROUTINE prepare_q()
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : bg, tpiba
  USE mp_pools, ONLY : npool, inter_pool_comm
  USE mp_world, ONLY : world_comm, nproc
  USE constants, ONLY : pi
  USE io_global, ONLY : stdout
  USE exx, ONLY : ecutfock
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE cell_base, ONLY : at
  USE control_ph, ONLY : current_iq
  USE noncollin_module, ONLY : npol
  USE uspp, ONLY : nkb, okvan
  USE exx, ONLY : dfftt
  USE mp, ONLY : mp_circular_shift_left, mp_barrier
  !
  USE sctk_val, ONLY : gindx, gq2, ngv, nqbz, nk_p, k0_p, nk_p_max, &
  &                    wfc, wfcq, becwfc, becwfcq, wscr, &
  &                    ngv0, ngv1, lsf, nmf, nbnd_p, nbnd_p_max
  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik_p, ik_g, jk_p, jk_g, ib, i1, i2, i3, ikv(3), jkv(3), g0(3), ir(3), ifft, ig, igv(3), &
  &          my_k0_p, my_nk_p, jkindx(nqbz), ipe, iqv(3), ipol
  REAL(dp) :: gv(3), theta, gq20, RAM_V, RAM_rho
  COMPLEX(dp) :: phase(dfftt%nnr)
  COMPLEX(DP),ALLOCATABLE :: wfctmp(:,:,:,:), becwfctmp(:,:,:,:)
  !
  ALLOCATE(wfctmp(dfftt%nnr,nbnd_p_max,npol,nk_p_max), becwfctmp(nkb,nbnd_p_max,npol,nk_p_max))
  !
  ! |G+q|^2
  !
  ifft = 0
  ngv = 0
  DO i3 = 1, dfftt%nr3
     DO i2 = 1, dfftt%nr2
        DO i1 = 1, dfftt%nr1
           !
           ifft = ifft + 1
           !
           igv(1:3) = (/i1, i2, i3/) - 1
           WHERE(igv(1:3)*2 >= (/dfftt%nr1, dfftt%nr2, dfftt%nr3/)) &
           &  igv(1:3) = igv(1:3) - (/dfftt%nr1, dfftt%nr2, dfftt%nr3/)
           gv(1:3) = MATMUL(bg(1:3,1:3), REAL(igv(1:3), dp)) * tpiba
           gv(1:3) = gv(1:3) - x_q(1:3, current_iq) * tpiba
           !
           gq20 = DOT_PRODUCT(gv(1:3), gv(1:3))
           !
           IF(ecutfock < 1e-10_dp .OR. gq20 < ecutfock) THEN
              !
              ngv = ngv + 1
              gq2(ngv) = gq20 / (8.0_dp * pi)
              gindx(ngv) = ifft
              !
           END IF
           !
        END DO ! i1 = 1, dfftt%nr1
     END DO ! i2 = 1, dfftt%nr2
  END DO ! i3 = 1, dfftt%nr3
  !
  WRITE(stdout,'(9x,"Number of PWs for W : ",i0)') ngv
  !
  CALL divide(world_comm, ngv, ngv0, ngv1)
  !
  ! RAM size estimation
  !
  RAM_V = REAL(ngv,dp)**2*REAL(nmf+1)
  RAM_rho = REAL(ngv + (ngv1-ngv0),dp)*REAL(nbnd_p_max,dp)*REAL(nk_p_max,dp)*REAL(nproc,dp)
  !
  IF(lsf>0) THEN
     IF(npol == 2) THEN
        RAM_V = RAM_V*4.0_dp
        RAM_rho = RAM_rho*4.0_dp
     ELSE
        RAM_V = RAM_V*2.0_dp
     END IF
  END IF
  !
  WRITE(stdout,'(9x,"Total RAM for Vscr : ",e10.2," GB")') RAM_V*16.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for rho  : ",e10.2," GB")') RAM_rho*16.0e-9_dp
  CALL mp_barrier(world_comm)
  IF(lsf>0) THEN
     ALLOCATE(wscr(ngv, ngv0:ngv1, 0:nmf, 2*npol))
  ELSE
     ALLOCATE(wscr(ngv, ngv0:ngv1, 0:nmf, 1))
  END IF
  !
  ! Prepare wave functions with phase shift
  !
  iqv(1:3) = NINT(MATMUL(x_q(1:3,current_iq), at(1:3,1:3)) * REAL((/nq1, nq2, nq3/), dp) - 0.5_dp)
  !
  DO ik_g = 1, nqbz
     !
     ikv(1) = (ik_g - 1) / (nq3*nq2)
     ikv(2) = (ik_g - 1 - ikv(1)*nq2*nq3) / nq3
     ikv(3) =  ik_g - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
     !
     WHERE(ikv(1:3)*2 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
     !
     jkv(1:3) = ikv(1:3) + iqv(1:3)
     jkv(1:3) = MODULO(jkv(1:3), (/nq1,nq2,nq3/))
     WHERE(jkv(1:3)*2 + 1 >= (/nq1,nq2,nq3/)) jkv(1:3) = jkv(1:3) - (/nq1,nq2,nq3/)
     !
     ! G0 = k' - (k + q)
     !
     g0(1:3) = (jkv(1:3) - ikv(1:3) - iqv(1:3)) / (/nq1,nq2,nq3/)
     !
     jkv(1:3) = MODULO(jkv(1:3), (/nq1,nq2,nq3/))
     jk_g = 1 + jkv(3) + jkv(2)*nq3 + jkv(1)*nq2*nq3
     jkindx(ik_g) = jk_g
     !
     IF(ik_g <= k0_p .OR. k0_p+nk_p < ik_g) CYCLE
     !
     ik_p = ik_g - k0_p
     !
     ! phi_{k}(r) = exp(i k r) u_k(r)
     ! phi_{k+G0}(r) = exp(i (k+G0) r) u_{k+G0}(r) = phi_{k}(r)
     ! u_{k+G0}(r) = exp(-i G0 r) u_k(r)
     !
     ! u_{k+q}(r) = u_{k'-G0}(r) = exp(i G0 r) u_{k'}(r)
     !
     ig = 0
     DO i3 = 1, dfftt%nr3
        DO i2 = 1, dfftt%nr2
           DO i1 = 1, dfftt%nr1
              !
              ig = ig + 1
              !
              ir(1:3) = (/i1, i2, i3/) - 1
              ir(1:3) = ir(1:3) * g0(1:3)
              !             
              theta = SUM(REAL(ir(1:3), dp) / REAL((/dfftt%nr1, dfftt%nr2, dfftt%nr3/), dp))
              theta = 2.0_dp * pi * theta
              !
              phase(ig) = CMPLX(COS(theta), SIN(theta), KIND=dp)
              !
           END DO ! i1
        END DO ! i2
     END DO ! i3
     !
     DO ipol = 1, npol
        DO ib = 1, nbnd_p
           !
           ! This phase should be applied as phi(2) * phase.
           ! But phi(2) will be rearranged below.
           ! Therefore, it is applied as
           ! phi(1)^* phi(2) phase = (phi(1) phase^*)^* phi(2)
           !
           wfcq( 1:dfftt%nnr,ib,ipol,ik_p, 1) = &
           & wfc(1:dfftt%nnr,ib,ipol,ik_p, 1) * CONJG(phase(1:dfftt%nnr))
           wfcq( 1:dfftt%nnr,ib,ipol,ik_p, 2) = &
           & wfc(1:dfftt%nnr,ib,ipol,ik_p, 2)
           !
        END DO ! ipol
     END DO ! ib
     !
  END DO ! ik
  !
  IF(okvan) becwfcq(1:nkb,1:nbnd_p,1:npol,1:nk_p, 1:2) &
  &        = becwfc(1:nkb,1:nbnd_p,1:npol,1:nk_p, 1:2)
  !
  wfctmp(1:dfftt%nnr,1:nbnd_p,1:npol,1:nk_p) = wfcq(1:dfftt%nnr,1:nbnd_p,1:npol,1:nk_p,2)
  wfcq(  1:dfftt%nnr,1:nbnd_p,1:npol,1:nk_p,2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  IF (okvan) THEN
     becwfctmp(1:nkb,1:nbnd_p,1:npol,1:nk_p) = becwfcq(1:nkb,   1:nbnd_p,1:npol,1:nk_p,2)
     becwfcq(  1:nkb,1:nbnd_p,1:npol,1:nk_p,2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  END IF
  !
  my_k0_p = k0_p
  my_nk_p = nk_p
  !
  DO ipe = 1, npool
     !
     CALL mp_circular_shift_left(nk_p, 1, inter_pool_comm )
     CALL mp_circular_shift_left(k0_p,  1, inter_pool_comm )
     CALL circular_shift_wrapper_c(    dfftt%nnr*npol*nbnd_p_max,nk_p_max,inter_pool_comm,wfctmp)
     IF(okvan) CALL circular_shift_wrapper_c(nkb*npol*nbnd_p_max,nk_p_max,inter_pool_comm,becwfctmp)
     !
     DO ik_p = 1, my_nk_p
        !
        ik_g = ik_p + my_k0_p
        jk_g = jkindx(ik_g)
        IF(jk_g <= k0_p .OR. k0_p+nk_p < jk_g) CYCLE
        !
        jk_p = jk_g - k0_p
        !
        wfcq(    1:dfftt%nnr,1:nbnd_p,1:npol,ik_p,2) = &
        & wfctmp(1:dfftt%nnr,1:nbnd_p,1:npol,jk_p)
        IF (okvan) THEN
           becwfcq(    1:nkb,1:nbnd_p,1:npol,ik_p,2) = &
           & becwfctmp(1:nkb,1:nbnd_p,1:npol,jk_p)
        END IF
        !
     END DO
     !
  END DO
  !
  DEALLOCATE(wfctmp, becwfctmp)
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
  USE mp_world, ONLY : world_comm
  USE mp_pools, ONLY : inter_pool_comm, npool
  USE mp, ONLY : mp_sum, mp_barrier, mp_circular_shift_left
  USE wvfct, ONLY : et, nbnd
  USE klist, ONLY : nks
  USE control_ph, ONLY : current_iq
  USE io_global, ONLY : stdout
  !
  USE sctk_tetra, ONLY : tetraweight, interpol_indx
  USE sctk_val, ONLY : nbnd_p, nbnd_p_max, nk_p, nk_p_max, k0_p, nmf
  !
  IMPLICIT NONE
  !
  COMPLEX(dp),INTENT(OUT) :: wght((nmf+1)*nbnd*nbnd_p_max,nk_p_max) !< integration weight
  !
  INTEGER :: ntetra0, ntetra1, nks0, it, ik, ii, dik(3), ikv(3), ik_p, &
  &          indx2(20, ntetra), indx3(20 * ntetra), kintp(8), ikq, iproc
  REAL(dp) :: kv(3), wintp(8), RAM_wghtd
  COMPLEX(dp),ALLOCATABLE :: wghtd(:,:)
  !
  CALL start_clock("fermi_factor")
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
  CALL divide(inter_pool_comm, ntetra,ntetra0,ntetra1)
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
  RAM_wghtd = REAL(nmf+1, dp)*REAL(nbnd, dp)*REAL(nbnd_p, dp)*REAL(nks0, dp)
  CALL mp_sum(RAM_wghtd, world_comm)
  WRITE(stdout,'(9x,"Total RAM for Weight on dense grid : ",e10.2," GB")') RAM_wghtd*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wghtd((nmf+1)*nbnd*nbnd_p,nks0))
  !
  CALL tetraweight(nks0,indx2,wghtd)
  !
  ! Interpolation of weight
  !
  wght(1:(nmf+1)*nbnd*nbnd_p_max,1:nk_p_max) = 0.0_dp
  !
  DO iproc = 1, npool
     !
     CALL mp_circular_shift_left(wght,1,inter_pool_comm)
     CALL mp_circular_shift_left(nk_p,1,inter_pool_comm)
     CALL mp_circular_shift_left(k0_p,1,inter_pool_comm)
     !
     DO ik = 1, nks0
        !
        ikv(1) = (indx3(ik) - 1) / (nk3*nk2)
        ikv(2) = (indx3(ik) - 1 - ikv(1)*nk2*nk3) / nk3
        ikv(3) =  indx3(ik) - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
        !
        kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
        CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
        !
        DO ii = 1, 8
           IF(k0_p < kintp(ii) .AND. kintp(ii) <= k0_p + nk_p) THEN
              !
              ik_p = kintp(ii) - k0_p
              !
              wght(     1:(nmf+1)*nbnd*nbnd_p, ik_p) &
              & = wght( 1:(nmf+1)*nbnd*nbnd_p, ik_p) &
              & + wghtd(1:(nmf+1)*nbnd*nbnd_p, ik) * wintp(ii)
              !
           END IF
        END DO ! ii
        !
     END DO ! ik
  END DO ! iproc
  !
  DEALLOCATE(wghtd)
  !
  CALL stop_clock("fermi_factor")
  !
END SUBROUTINE fermi_factor
!>
!> Compute screened interaction
!>
SUBROUTINE make_scrn()
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE mp_world, ONLY : nproc
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE us_exx, ONLY : addusxx_r
  USE mp, ONLY : mp_circular_shift_left, mp_sum, mp_barrier
  USE mp_world, ONLY : world_comm
  USE io_global, ONLY : stdout
  !
  USE sctk_val, ONLY : ngv, nmf, nbnd_p_max, nbnd_p, nk_p, &
  &                    wscr, lsf, ngv0, ngv1, nk_p_max
  !
  IMPLICIT NONE
  !
  INTEGER :: ik_p, ibnd_p, jbnd_g, ikib, imf, ipe, ipol, npol2, nwscr
  REAL(DP) :: RAM_wght
  !
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:,:), wght(:,:,:,:)
  !
  CALL start_clock("make_scrn")
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
  wscr(1:ngv, ngv0:ngv1, 0:nmf, 1:nwscr) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Calc f * (1 - f') / (e - e' + iw)
  !
  RAM_wght = REAL(nmf+1,dp)*REAL(nbnd,DP)*REAL(nbnd_p_max)*REAL(nk_p_max)
  CALL mp_sum(RAM_wght, world_comm)
  WRITE(stdout,'(9x,"Total RAM for weight on coarse grid : ",e10.2," GB")') RAM_wght*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wght(0:nmf,nbnd,nbnd_p_max,nk_p_max))
  CALL fermi_factor(wght)
  !
  ! Average the weights for degenerated states
  !
  CALL average_wght_kel(nmf+1, wght(0:nmf,1:nbnd,1:nbnd_p_max,1:nk_p_max))
  !
  ! Calc. Chi
  !
  ALLOCATE(rho1(ngv,nbnd_p_max*nk_p_max,npol2), rho2(ngv0:ngv1,nbnd_p_max*nk_p_max,npol2))
  !
  DO jbnd_g = 1, nbnd
     !
     CALL calc_rhog(npol2, jbnd_g, rho1)
     !
     DO ipe = 1, nproc
        !
        CALL mp_circular_shift_left(nbnd_p, 1, world_comm )
        CALL mp_circular_shift_left(nk_p, 1, world_comm )
        CALL circular_shift_wrapper_c(ngv,nbnd_p_max*nk_p_max*npol2,world_comm,rho1)
        CALL circular_shift_wrapper_c(nmf+1,nbnd_p_max*nbnd*nk_p_max,world_comm,wght)
        !
        DO imf = 0, nmf
           !
           ikib = 0
           DO ik_p = 1, nk_p
              DO ibnd_p = 1, nbnd_p
                 ikib = ikib + 1
                 rho2(ngv0:ngv1,ikib,1:npol2) = REAL(wght(imf,jbnd_g,ibnd_p,ik_p),dp) &
                 &                            * CONJG(rho1(ngv0:ngv1,ikib,1:npol2))
              END DO
           END DO
           !
           CALL start_clock("zgemm(make_scrn)")
           !
           DO ipol = 1, npol2
              CALL zgemm("N", "T", ngv, ngv1-ngv0+1, ikib, &
              &  (1.0_dp, 0.0_dp), rho1(1:ngv,            1:ikib, ipol), ngv, &
              &                    rho2(       ngv0:ngv1, 1:ikib, ipol), ngv1-ngv0+1, &
              &  (1.0_dp, 0.0_dp), wscr(1:ngv, ngv0:ngv1,    imf, ipol), ngv  )
           END DO
           !
           CALL stop_clock("zgemm(make_scrn)")
           !
        END DO ! imf = 0, nmf
        !
     END DO ! ipe = 1, nproc
     !
  END DO
  !
  IF(nwscr == 2) &
  &  wscr(1:ngv, ngv0:ngv1, 0:nmf, 2) = wscr(1:ngv, ngv0:ngv1, 0:nmf, 1)
  !
  DEALLOCATE(rho1, rho2, wght)
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
  USE mp_world, ONLY : world_comm, nproc
  USE mp, ONLY : mp_sum, mp_circular_shift_left, mp_barrier
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE us_exx, ONLY : addusxx_r
  USE io_global, ONLY : stdout
  !
  USE sctk_val, ONLY : gq2, ngv, nmf, nbnd_p, nbnd_p_max, &
  &                    ngv0, ngv1, Kel, Wscr, lsf, nk_p, nk_p_max
  !
  IMPLICIT NONE
  !
  INTEGER :: ik_p, ibnd_p, jbnd_g, imf, ipe, ipol, jpol, npol2, nwscr, &
  &          ikib, nKel
  REAL(dp) :: RAM_Kel
  !
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:), Kel0(:,:), Kel_temp(:,:,:,:)
  COMPLEX(DP),EXTERNAL :: zdotc
  !
  CALL start_clock("make_kel")
  !
  IF(lsf>0) THEN
     IF(npol == 2) THEN
        npol2 = 4
        nwscr = 4
     ELSE
        npol2 = 1
        nwscr = 2
     END IF
     nKel = 2
  ELSE
     npol2 = 1
     nwscr = 1
     nKel = 1
  END IF
  !
  RAM_Kel = REAL(nmf+2,dp)*REAL(nbnd,DP)*REAL(nbnd_p_max)*REAL(nk_p_max)*REAL(nKel,DP)
  CALL mp_sum(RAM_Kel, world_comm)
  WRITE(stdout,'(9x,"Total RAM for Kel : ",e10.2," GB")') RAM_Kel*8.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(Kel(0:nmf+1,nbnd,nbnd_p_max,nk_p_max,nKel))
  Kel(0:nmf + 1,1:nbnd,1:nbnd_p_max,1:nk_p_max,1:nKel) = 0.0_dp
  !
  ALLOCATE(rho1(ngv,nbnd_p_max*nk_p_max,npol2), rho2(nbnd_p_max*nk_p_max,ngv0:ngv1), &
  &        Kel0(nbnd_p_max*nk_p_max,nwscr))
  !
  DO jbnd_g = 1, nbnd
     !
     CALL calc_rhog(npol2, jbnd_g, rho1)
     !
     ! Infinite frequency -> Bare Coulomb
     !
     ikib = 0
     DO ik_p = 1, nk_p
        DO ibnd_p = 1, nbnd_p
           ikib = ikib + 1
           Kel(nmf + 1,jbnd_g,ibnd_p,ik_p,1) = SUM(REAL(rho1(1:ngv,ikib,1) &
           &                                    * CONJG(rho1(1:ngv,ikib,1)), dp) &
           &                                    /        gq2(1:ngv))  / omega
        END DO ! ibnd_p
     END DO ! ik_p
     !
     DO ipe = 1, nproc
        !
        CALL mp_circular_shift_left(nk_p, 1, world_comm )
        CALL mp_circular_shift_left(nbnd_p,  1, world_comm )
        CALL circular_shift_wrapper_c(ngv,nk_p_max*nbnd_p_max*npol2,world_comm,rho1)
        IF(lsf>0) THEN
           CALL circular_shift_wrapper_r(nmf + 2,nbnd*nbnd_p_max*nk_p_max*2,world_comm,Kel)
        ELSE
           CALL circular_shift_wrapper_r(nmf + 2,nbnd*nbnd_p_max*nk_p_max,world_comm,Kel)
        END IF
        !
        DO imf = 0, nmf
           !
           DO ipol = 1, nwscr
              !
              IF(nwscr == 2) THEN
                 jpol = 1
              ELSE
                 jpol = ipol
              END IF
              !
              CALL start_clock("zgemm(make_kel)")
              !
              CALL zgemm("T", "N", nbnd_p*nk_p, ngv1-ngv0+1, ngv, &
              &          (1.0_dp, 0.0_dp), rho1(1:ngv,1:nbnd_p*nk_p,jpol), ngv, &
              &                            Wscr(1:ngv,                      ngv0:ngv1, imf,ipol),ngv, &
              &          (0.0_dp, 0.0_dp), rho2(      1:nbnd_p_max*nk_p_max,ngv0:ngv1), nbnd_p_max*nk_p_max)
              !
              CALL stop_clock("zgemm(make_kel)")
              !
              DO ikib = 1, nbnd_p*nk_p
                 Kel0(ikib,ipol) = zdotc(ngv1-ngv0+1, &
                 &                       rho2(ikib,ngv0), nbnd_p_max*nk_p_max, &
                 &                       rho1(ngv0:ngv1,ikib,jpol), 1)
              END DO !ikib
              !
           END DO ! ipol = 1, nwscr
           !
           ikib = 0
           DO ik_p = 1, nk_p
              DO ibnd_p = 1, nbnd_p
                 ikib = ikib + 1
                 !
                 Kel(imf,   jbnd_g,ibnd_p,ik_p,1) = Kel(imf,jbnd_g,ibnd_p,ik_p,1) + REAL(Kel0(ikib,1), dp) / omega
                 !
                 IF(nwscr == 2) THEN
                    Kel(imf,jbnd_g,ibnd_p,ik_p,2) = Kel(imf,jbnd_g,ibnd_p,ik_p,2) + REAL(Kel0(ikib,2), dp) / omega
                 ELSE IF(nwscr == 4) THEN
                    Kel(imf,jbnd_g,ibnd_p,ik_p,2) = Kel(imf,jbnd_g,ibnd_p,ik_p,2) + REAL(SUM(Kel0(ikib,2:4)), dp) / omega
                 END IF
                 !
              END DO ! ibnd_p
           END DO ! ik_p
           !
        END DO ! imf = 0, nmf
        !
     END DO ! ipe = 1, nproc
     !
  END DO ! jbnd
  !
  DEALLOCATE(rho1, rho2, Kel0)
  !
  ! Average the weights for degenerated states
  !
  ALLOCATE(Kel_temp(0:nmf + 1,nbnd,nbnd_p_max,nk_p_max))
  IF(lsf>0) THEN
     Kel_temp(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p) = CMPLX(Kel(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p,1), &
     &                                                  Kel(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p,2), DP)
  ELSE
     Kel_temp(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p) = CMPLX(Kel(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p,1), 0.0_dp, DP)
  END IF
  CALL average_wght_kel(nmf+2, Kel_temp(0:nmf+1,1:nbnd,1:nbnd_p_max,1:nk_p_max))
  Kel(          0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p,1) = REAL( Kel_temp(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p), DP)
  IF(lsf>0) Kel(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p,2) = AIMAG(Kel_temp(0:nmf + 1,1:nbnd,1:nbnd_p,1:nk_p))
  !
  DEALLOCATE(Kel_temp)
  !
  CALL stop_clock("make_kel")
  !
END SUBROUTINE make_Kel
!>
!>
!>
SUBROUTINE average_wght_kel(nmf0,wght)
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE mp_pools, ONLY : intra_pool_comm, nproc_pool
  USE mp, ONLY : mp_circular_shift_left
  USE mp, ONLY : mp_barrier !debug
  USE mp_world, ONLY : world_comm ! debug
  !
  USE sctk_val, ONLY : nbnd_p, nbnd_p_max, bnd0_p, nk_p, nk_p_max, k0_p, degen, ndegen
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nmf0
  COMPLEX(DP),INTENT(INOUT) :: wght(nmf0,nbnd,nbnd_p_max,nk_p_max)
  !
  INTEGER :: ik_p, ik_g, ibnd, my_nbnd_p, my_bnd0_p, iproc, idegen
  COMPLEX(DP) :: wght_ave(nmf0,nbnd_p_max)
  COMPLEX(DP),ALLOCATABLE :: wght_tmp(:,:,:,:)
  !
  ALLOCATE(wght_tmp(nmf0,nbnd_p_max,nbnd,nk_p_max))
  !
  DO ik_p = 1, nk_p
     !
     ik_g = ik_p + k0_p
     !
     DO idegen = 1, ndegen(ik_g,2)
        wght_ave(1:nmf0,1:nbnd_p) = SUM(wght(1:nmf0,degen(1,idegen,ik_g,2): &
        &                                           degen(2,idegen,ik_g,2),1:nbnd_p,ik_p), 2) &
        &                      / REAL(degen(2,idegen,ik_g,2) - degen(1,idegen,ik_g,2) + 1, DP)
        DO ibnd = degen(1,idegen,ik_g,2), degen(2,idegen,ik_g,2)
           wght(1:nmf0,ibnd,1:nbnd_p,ik_p)  = wght_ave(1:nmf0,1:nbnd_p)
        END DO
     END DO
     !
  END DO
  !
  wght_tmp(1:nmf0,1:nbnd_p_max,1:nbnd,1:nk_p_max) = 0.0_dp
  my_nbnd_p = nbnd_p
  my_bnd0_p = bnd0_p
  DO iproc = 1, nproc_pool
     !
     CALL mp_circular_shift_left(bnd0_p, 1, intra_pool_comm )
     CALL mp_circular_shift_left(nbnd_p, 1, intra_pool_comm )
     CALL circular_shift_wrapper_c(nmf0,nbnd*nbnd_p_max*nk_p_max,intra_pool_comm,wght)
     !
     wght_tmp(1:nmf0,          1:          my_nbnd_p,bnd0_p+1:bnd0_p+nbnd_p,1:nk_p_max) &
     & = wght(1:nmf0,my_bnd0_p+1:my_bnd0_p+my_nbnd_p,       1:       nbnd_p,1:nk_p_max)
     !
  END DO
  !
  DO ik_p = 1, nk_p
     !
     ik_g = ik_p + k0_p
     !
     DO idegen = 1, ndegen(ik_g,1)
        wght_ave(1:nmf0,1:my_nbnd_p) = SUM(wght_tmp(1:nmf0,1:my_nbnd_p,degen(1,idegen,ik_g,1): &
        &                                                            degen(2,idegen,ik_g,1),ik_p), 3) &
        &                      / REAL(degen(2,idegen,ik_g,1) - degen(1,idegen,ik_g,1) + 1, DP)
        DO ibnd = degen(1,idegen,ik_g,1), degen(2,idegen,ik_g,1)
           wght_tmp(1:nmf0,1:my_nbnd_p,ibnd,ik_p)  = wght_ave(1:nmf0,1:my_nbnd_p)
        END DO
     END DO
     !
  END DO
  !
  DO iproc = 1, nproc_pool
     !
     CALL mp_circular_shift_left(bnd0_p, 1, intra_pool_comm )
     CALL mp_circular_shift_left(nbnd_p, 1, intra_pool_comm )
     CALL circular_shift_wrapper_c(nmf0,nbnd*nbnd_p_max*nk_p_max,intra_pool_comm,wght)
     !
     wght(        1:nmf0,my_bnd0_p+1:my_bnd0_p+my_nbnd_p,       1:       nbnd_p,1:nk_p_max) &
     & = wght_tmp(1:nmf0,          1:          my_nbnd_p,bnd0_p+1:bnd0_p+nbnd_p,1:nk_p_max)
     !
  END DO
  !
  DEALLOCATE(wght_tmp)
  !
  CALL mp_barrier(world_comm)
  !
END SUBROUTINE average_wght_kel
!>
!>
!>
SUBROUTINE calc_rhog(npol2, jbnd_g, rho1)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE exx, ONLY : dfftt
  USE us_exx, ONLY : addusxx_r
  USE uspp, ONLY : nkb, okvan
  USE mp_pools, ONLY : intra_pool_comm
  USE mp, ONLY : mp_max, mp_circular_shift_left, mp_sum
  !
  USE sctk_val, ONLY : gindx, ngv, wfcq, becwfcq, nk_p, nk_p_max, nbnd_p, nbnd_p_max, bnd0_p
  !
  IMPLICIT none
  !
  INTEGER,INTENT(IN) :: npol2, jbnd_g
  COMPLEX(DP),INTENT(OUT) :: rho1(ngv,nbnd_p_max*nk_p_max,npol2)
  !
  INTEGER:: ibnd_p, jbnd_p, ik_p, ikib, ipol
  COMPLEX(dp) :: rho0(dfftt%nnr,4), rho4(dfftt%nnr)
  COMPLEX(dp),ALLOCATABLE :: wfc2(:,:,:), becwfc2(:,:,:)
  !
  ALLOCATE(wfc2(dfftt%nnr,npol,nk_p_max))
  IF(okvan) ALLOCATE(becwfc2(nkb,npol,nk_p_max))
  !
  IF(bnd0_p < jbnd_g .AND. jbnd_g <= bnd0_p+nbnd_p) THEN
     jbnd_p = jbnd_g - bnd0_p
     wfc2(1:dfftt%nnr,1:npol,1:nk_p_max) = wfcq(1:dfftt%nnr,jbnd_p,1:npol,1:nk_p_max,2)
     IF(okvan) becwfc2(1:nkb,1:npol,1:nk_p_max) = becwfcq(1:nkb,jbnd_p,1:npol,1:nk_p_max,2)
  ELSE
     wfc2(1:dfftt%nnr,1:npol,1:nk_p_max) = 0.0_dp
     IF(okvan) becwfc2(1:nkb,1:npol,1:nk_p_max) = 0.0_dp
  END IF
  CALL mp_sum(wfc2, intra_pool_comm)
  IF(okvan) CALL mp_sum(becwfc2, intra_pool_comm)
  !
  ikib = 0
  DO ik_p = 1, nk_p
     !
     DO ibnd_p = 1, nbnd_p
        !
        ikib = ikib + 1
        !
        ! debug
        !rho4(1:dfftt%nnr) = conjg(wfcq(1:dfftt%nnr,ibnd_p,1,ik_p, 1)) * wfcq(1:dfftt%nnr,ibnd_p,1,ik_p, 1)
        !write(*,'(a,2f15.8)',advance="no") "debug1", sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
        !IF(okvan) then
        !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) / omega
        !   CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,1,ik_p,1), &
        !   &                                 becwfcq(1:nkb,ibnd_p,1,ik_p,1))
        !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) * omega
        !   write(*,'(2f15.8)',advance="no") sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
        !END IF
        !rho4(1:dfftt%nnr) = conjg(wfcq(1:dfftt%nnr,ibnd_p,1,ik_p, 2)) * wfcq(1:dfftt%nnr,ibnd_p,1,ik_p, 2)
        !write(*,'(a,2f15.8)',advance="no") "debug1", sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
        !IF(okvan) then
        !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) / omega
        !   CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,1,ik_p,2), &
        !   &                                 becwfcq(1:nkb,ibnd_p,1,ik_p,2))
        !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) * omega
        !   write(*,'(2f15.8)') sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
        !END IF
        ! end debug
        !
        CALL start_clock("wfc*wfc=rho")
        !
        rho0(1:dfftt%nnr, 1:npol2) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        DO ipol = 1, npol
           !
           rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr,ibnd_p,ipol,ik_p,1)) &
           &                       * wfc2(1:dfftt%nnr,       ipol,ik_p) / omega
           IF(okvan) CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,ipol,ik_p,1), &
           &                                           becwfc2(1:nkb,       ipol,ik_p))
           rho0(1:dfftt%nnr,1) = rho0(1:dfftt%nnr,1) + rho4(1:dfftt%nnr) * omega
           !
           IF(npol2 == 4) THEN
              !
              rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr,ibnd_p,  ipol,ik_p,1)) &
              &                       * wfc2(1:dfftt%nnr,       3-ipol,ik_p) / omega
              IF(okvan) CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,  ipol,ik_p,1), &
              &                                           becwfc2(1:nkb,       3-ipol,ik_p))
              rho0(1:dfftt%nnr,2) = rho0(1:dfftt%nnr,2) + rho4(1:dfftt%nnr) * omega
              !
              rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr,ibnd_p,  ipol,ik_p,1)) &
              &                       * wfc2(1:dfftt%nnr,       3-ipol,ik_p) / omega
              IF(okvan) CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,  ipol,ik_p,1), &
              &                                           becwfc2(1:nkb,       3-ipol,ik_p))
              rho0(1:dfftt%nnr,3) = rho0(1:dfftt%nnr,3) &
              &                   + rho4(1:dfftt%nnr) * REAL(3-2*ipol, dp) * omega
              !
              rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr,ibnd_p,ipol,ik_p,1)) &
              &                       * wfc2(1:dfftt%nnr,       ipol,ik_p) / omega
              IF(okvan) CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd_p,ipol,ik_p,1), &
              &                                           becwfc2(1:nkb,       ipol,ik_p))
              rho0(1:dfftt%nnr,4) = rho0(1:dfftt%nnr,4) &
              &                   + rho4(1:dfftt%nnr) * REAL(3-2*ipol, dp) * omega
              !
           END IF
           !
        END DO ! ipol
        !
        CALL stop_clock("wfc*wfc=rho")
        !
        CALL start_clock("fft[rho]")
        !
        DO ipol = 1, npol2
           CALL cfft3d (rho0(1:dfftt%nnr,ipol), &
           & dfftt%nr1, dfftt%nr2, dfftt%nr3, dfftt%nr1, dfftt%nr2, dfftt%nr3, 1, -1)
           !
           rho1(1:ngv, ikib,ipol) = rho0(gindx(1:ngv),ipol)
        END DO ! ipol
        !
        CALL stop_clock("fft[rho]")
        !
     END DO ! ibnd_p
  END DO ! ik_p
  !
END SUBROUTINE calc_rhog
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
  USE mp, ONLY : mp_sum, mp_circular_shift_left
  !
  USE sctk_val, ONLY : nmf, mf, Kel, lsf, bnd0_p, nbnd_p, &
  &                    k0_p, nk_p
  !
  IMPLICIT NONE
  !
  INTEGER :: ntetra0, ntetra1, ii, ikv(3), ik, ik_p, imf, kintp(8), nkel, ikel, &
  &          b_low_p, b_high_p, b_low_g, b_high_g
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
  !
  b_low_g = MAX(b_low, bnd0_p+1)
  b_high_g = MIN(b_high, bnd0_p+nbnd_p)
  b_low_p = b_low_g - bnd0_p
  b_high_p = b_high_g - bnd0_p
  !
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
              !
              IF(k0_p < kintp(ii) .AND. kintp(ii) <= k0_p + nk_p) THEN
                 !
                 ik_p = kintp(ii) - k0_p
                 !
                 mu(imf,ikel) = mu(imf,ikel) + &
                 & SUM(Kel(imf,b_low:b_high,b_low_p:b_high_p,       ik_p, ikel) &
                 &     * wght( b_low:b_high,b_low_g:b_high_g,ik)) * wintp(ii)
                 !
              END IF
              !
           END DO ! ii = 1, 8
        END DO ! imf = 0, nmf + 1
     END DO ! ikel = 1, nkel
     !
  END DO ! ik = 1, nks
  !
  CALL mp_sum(mu, world_comm)
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
  USE sctk_val, ONLY : Kel, nqbz, nmf, lsf, mf, nci, nk_p, nbnd_p
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jb, imf, jmf, ikel, nkel
  REAL(dp) :: coef(nci,nci), Kel0(nci), Kel1(0:nmf+1), cheb(1:nci,0:nmf+1), &
  &           err(0:nmf+1), aveer(0:nmf+1), maxer(0:nmf+1), mf0, x0
  !
  IF(lsf>0) THEN
     nkel = 2
  ELSE
     nkel = 1
  END IF
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
     DO ik = 1, nk_p
        DO ib = 1, nbnd_p
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
END SUBROUTINE chebyshev_interpol
!>
!> Output to file
!>
SUBROUTINE write_Kel()
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE mp_world, ONLY : mpime, world_comm, nproc
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE cell_base, ONLY : at
  USE control_ph, ONLY : current_iq
  USE mp, ONLY : mp_circular_shift_left
  USE io_files, ONLY : prefix, tmp_dir
  !
  USE sctk_val, ONLY : Kel, nqbz, lsf, nci, bnd0_p, nbnd_p, nbnd_p_max, &
  &                    k0_p, nk_p, nk_p_max, nmf
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, nkel, iproc
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  REAL(DP),ALLOCATABLE :: Kel_all(:,:,:,:,:)
  !
  IF(lsf>0) THEN
     nkel = 2
  ELSE
     nkel = 1
  END IF
  IF(mpime==0) ALLOCATE(Kel_all(0:nci-1,nbnd,nbnd,nqbz,nkel))
  !
  DO iproc = 1, nproc
     CALL mp_circular_shift_left(nbnd_p, 1, world_comm )
     CALL mp_circular_shift_left(bnd0_p, 1, world_comm )
     CALL mp_circular_shift_left(k0_p, 1, world_comm )
     CALL mp_circular_shift_left(nk_p, 1, world_comm )
     CALL circular_shift_wrapper_r(nmf + 2,nbnd*nbnd_p_max*nk_p_max*nkel,world_comm,Kel)
     IF(mpime==0) Kel_all(0:nci-1,1:nbnd,bnd0_p+1:bnd0_p+nbnd_p,k0_p+1:k0_p+nk_p,1:nkel) &
     &              = Kel(0:nci-1,1:nbnd,       1:       nbnd_p,     1:     nk_p,1:nkel)
  END DO
  !
  IF(mpime == 0) THEN
     !
     OPEN(fo, file = TRIM(tmp_dir) // TRIM(prefix) // ".vel"//TRIM(int_to_char(current_iq)), &
     &    form = "unformatted")
     !
     WRITE(fo) nq1, nq2, nq3
     WRITE(fo) nbnd
     WRITE(fo) MATMUL(x_q(1:3, current_iq), at(1:3, 1:3))
     WRITE(fo) nci
     !
     WRITE(fo) Kel_all(0:nci-1,1:nbnd,1:nbnd,1:nqbz,1)
     IF(lsf>0) THEN
        WRITE(fo) Kel_all(0:nci-1,1:nbnd,1:nbnd,1:nqbz,2)
     END IF
     !
     CLOSE(fo)
     !
  END IF
  !
  IF(mpime==0) DEALLOCATE(Kel_all)
  DEALLOCATE(Kel)
  !
END SUBROUTINE write_Kel
!
END MODULE sctk_coulomb
