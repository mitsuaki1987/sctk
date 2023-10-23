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
  USE sctk_val, ONLY : gindx, gq2, nmf, nci, mf, kplusq, nqbz
  USE exx, ONLY : dfftt
  !
  IMPLICIT NONE
  !
  INTEGER :: imf
  REAL(dp) :: x0
  !
  nmf = 2*nci - 3
  !
  ALLOCATE(gq2(dfftt%nnr), gindx(dfftt%nnr), mf(nmf), kplusq(nqbz))
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
SUBROUTINE circular_shift_wrapper_c_nb(n1, n2, comm, wfc_s, wfc_r, req)
   !
   USE kinds, ONLY : DP
   USE mp, ONLY : mp_circular_shift_left_start
   !
   IMPLICIT NONE
   !
   INTEGER,INTENT(IN) :: n1, n2, comm
   COMPLEX(DP),INTENT(IN) :: wfc_s(n1,n2)
   COMPLEX(DP),INTENT(OUT) :: wfc_r(n1,n2)
   INTEGER,INTENT(INOUT) :: req(2)
   !
   CALL mp_circular_shift_left_start(wfc_s, wfc_r, 1, comm, req)
   !
 END SUBROUTINE circular_shift_wrapper_c_nb
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
  USE sctk_val, ONLY : gindx, gq2, ngv, nqbz, wfc, wfcq, becwfc, becwfcq, wscr, &
  &                    ngv0, ngv1, lsf, nmf, nb, nb_max, kplusq
  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, i1, i2, i3, ikv(3), jkv(3), g0(3), ir(3), ifft, ig, igv(3), &
  &          iqv(3), ipol
  REAL(dp) :: gv(3), theta, gq20, RAM_V, RAM_rho
  COMPLEX(dp) :: phase(dfftt%nnr)
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
  RAM_rho = REAL(dfftt%nnr + ngv + (ngv1-ngv0),dp)*REAL(nb_max,dp)*REAL(nqbz,dp)*REAL(nproc,dp)
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
     ! G0 = k' - (k + q)
     !
     g0(1:3) = (jkv(1:3) - ikv(1:3) - iqv(1:3)) / (/nq1,nq2,nq3/)
     !
     jkv(1:3) = MODULO(jkv(1:3), (/nq1,nq2,nq3/))
     kplusq(ik) = 1 + jkv(3) + jkv(2)*nq3 + jkv(1)*nq2*nq3
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
        !
        ! This phase should be applied as phi(2) * phase.
        ! But phi(2) will be rearranged below.
        ! Therefore, it is applied as
        ! phi(1)^* phi(2) phase = (phi(1) phase^*)^* phi(2)
        !
        DO ib = 1, nb(1)
           !
           wfcq( 1:dfftt%nnr, ib, ipol, ik, 1) = &
           & wfc(1:dfftt%nnr, ib, ipol, ik, 1) * CONJG(phase(1:dfftt%nnr))
           !
        END DO ! ib
        !
     END DO ! ipol
     !
     wfcq( 1:dfftt%nnr, 1:nb(2), 1:npol,        ik,  2) = &
     & wfc(1:dfftt%nnr, 1:nb(2), 1:npol, kplusq(ik), 2)
     !
     IF(okvan) THEN
        becwfcq(   1:nkb, 1:nb(1), 1:npol,        ik,  1) &
        & = becwfc(1:nkb, 1:nb(1), 1:npol,        ik,  1)
        becwfcq(   1:nkb, 1:nb(2), 1:npol,        ik,  2) &
        & = becwfc(1:nkb, 1:nb(2), 1:npol, kplusq(ik), 2)
     END IF
     !
  END DO ! ik
  !
END SUBROUTINE prepare_q
!>
!> Calculation of weight function 
!> @f$\frac{f(1-f')}{\varepsilon - \evarepsilon'}@f$
!>
SUBROUTINE fermi_factor(wght)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at
  USE disp,  ONLY : nq1, nq2, nq3, x_q
  USE start_k, ONLY : nk1, nk2, nk3
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum, mp_barrier
  USE wvfct, ONLY : et, nbnd
  USE klist, ONLY : nks
  USE control_ph, ONLY : current_iq
  USE io_global, ONLY : stdout
  !
  USE sctk_tetra, ONLY : tetraweight, interpol_indx
  USE sctk_val, ONLY : nb, nb_max, nqbz, nmf
  !
  IMPLICIT NONE
  !
  COMPLEX(dp),INTENT(OUT) :: wght(nmf+1, nqbz, nb_max, nb_max) !< integration weight
  !
  INTEGER :: ik, ii, dik(3), ikv(3), kintp(8), ikq
  REAL(dp) :: kv(3), wintp(8), RAM_wghtd
  COMPLEX(dp),ALLOCATABLE :: wghtd(:,:,:,:)
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
  RAM_wghtd = REAL(nmf+1, dp)*REAL(nb(1)*nb(2), dp)*REAL(nks, dp)
  CALL mp_sum(RAM_wghtd, world_comm)
  WRITE(stdout,'(9x,"Total RAM for Weight on dense grid : ",e10.2," GB")') RAM_wghtd*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wghtd(nmf+1, nb(2), nb(1), nks))
  !
  CALL tetraweight(wghtd)
  !
  ! Interpolation of weight
  !
  wght(1:nmf+1, 1:nqbz, 1:nb_max, 1:nb_max) = 0.0_dp
  !
  DO ik = 1, nks
     !
     ikv(1) = (ik - 1) / (nk3*nk2)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     !
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/), kv, kintp, wintp)
     !
     DO ii = 1, 8
        !
        wght(     1:nmf+1, kintp(ii), 1:nb(2), 1:nb(1)) &
        & = wght( 1:nmf+1, kintp(ii), 1:nb(2), 1:nb(1)) &
        & + wghtd(1:nmf+1,            1:nb(2), 1:nb(1), ik) * wintp(ii)
        !
     END DO ! ii
     !
  END DO ! ik
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
  USE kinds, ONLY : DP
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE exx, ONLY : dfftt
  USE mp, ONLY : mp_circular_shift_left, mp_sum, mp_barrier, mp_waitall
  USE mp_world, ONLY : world_comm
  USE io_global, ONLY : stdout
  USE uspp, ONLY : nkb, okvan
  !
  USE sctk_val, ONLY : ngv, nmf, nqbz, nb_max, wfcq, becwfcq, &
  &                    wscr, lsf, ngv0, ngv1, nb
  USE mp_kel, ONLY : nproc_band1, nproc_band2, band1_comm, band2_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, jbnd, imf, iproc1, iproc2, ipol, npol2, nwscr, jbnd_ik, req(6)
  REAL(DP) :: RAM_wght
  !
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:,:), wght(:,:,:,:), &
  &                          wfcq_r(:,:,:,:), becwfcq_r(:,:,:,:), wght_r(:,:,:,:)
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
  RAM_wght = REAL(nmf+1,dp)*REAL(nb_max**2,DP)*REAL(nqbz)
  CALL mp_sum(RAM_wght, world_comm)
  WRITE(stdout,'(9x,"Total RAM for weight on coarse grid : ",e10.2," GB")') RAM_wght*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wght(0:nmf,nqbz,nb_max,nb_max))
  CALL fermi_factor(wght)
  !
  ! Average the weights for degenerated states
  !
  CALL average_wght_kel(nmf+1, wght(0:nmf,1:nqbz,1:nb_max,1:nb_max))
  !
  ! Calc. Chi
  !
  ALLOCATE(rho1(ngv,nqbz*nb_max,npol2), rho2(ngv0:ngv1,nqbz*nb_max,npol2), &
  &        wfcq_r(dfftt%nnr, nb_max, npol, nqbz), wght_r(0:nmf,nqbz,nb_max,nb_max))
  IF(okvan) THEN
     ALLOCATE(becwfcq_r(nkb,nb_max, npol, nqbz))
  ELSE
     ALLOCATE(becwfcq_r(1,1, 1, 1))
  END IF
  !
  DO iproc1 = 1, nproc_band1
     !
     CALL mp_circular_shift_left(nb(1), 1, band1_comm )
     CALL circular_shift_wrapper_c(dfftt%nnr, nb_max*npol*nqbz, band1_comm, &
     &                             wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 1))
     IF(okvan) CALL circular_shift_wrapper_c(nkb, nb_max*npol*nqbz, band1_comm, &
     &                                       becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 1))
     CALL circular_shift_wrapper_c(nmf+1, nqbz*nb_max*nb_max, band1_comm, &
     &                             wght(0:nmf, 1:nqbz, 1:nb_max, 1:nb_max))
     !
     DO iproc2 = 1, nproc_band2
        !
        ! Non-blocking circular shift : They are used in the next loop.
        !
        CALL circular_shift_wrapper_c_nb(dfftt%nnr, nb_max*npol*nqbz, band2_comm, &
        &                             wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 2), &
        &                           wfcq_r(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz), req(1:2))
        IF(okvan) CALL circular_shift_wrapper_c_nb(nkb, nb_max*npol*nqbz, band2_comm, &
        &                                          becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 2), &
        &                                        becwfcq_r(1:nkb, 1:nb_max, 1:npol, 1:nqbz), req(5:6))
        CALL circular_shift_wrapper_c_nb(nmf+1, nqbz*nb_max*nb_max, band2_comm, &
        &                                wght(0:nmf, 1:nqbz, 1:nb_max, 1:nb_max), &
        &                              wght_r(0:nmf, 1:nqbz, 1:nb_max, 1:nb_max), req(3:4))
        !
        DO ibnd = 1, nb(1)
           !
           CALL calc_rhog(npol2, ibnd, rho1)
           !
           DO imf = 0, nmf
              !
              !$OMP PARALLEL DEFAULT(NONE) &
              !$OMP & SHARED(nb,nqbz,ngv0,ngv1,npol2,rho1,rho2,wght,imf,ibnd) &
              !$OMP & PRIVATE(jbnd_ik,jbnd,ik)
              !$OMP DO
              DO jbnd_ik = 1, nb(2)*nqbz
                 !
                 jbnd = 1 + (jbnd_ik-1) / nqbz
                 ik   = 1 +  jbnd_ik-1 -  nqbz*(jbnd-1)
                 !
                 rho2(ngv0:ngv1,jbnd_ik,1:npol2) = REAL(wght(imf,ik,jbnd,ibnd),dp) &
                 &                        * CONJG(rho1(ngv0:ngv1,jbnd_ik,1:npol2))
                 !
              END DO
              !$OMP END DO
              !$OMP END PARALLEL
              !
              CALL start_clock("zgemm(make_scrn)")
              !
              DO ipol = 1, npol2
                 CALL zgemm("N", "T", ngv, ngv1-ngv0+1, nqbz*nb(2), &
                 &  (1.0_dp, 0.0_dp), rho1(1:ngv,            1:nqbz*nb_max, ipol), ngv, &
                 &                    rho2(       ngv0:ngv1, 1:nqbz*nb_max, ipol), ngv1-ngv0+1, &
                 &  (1.0_dp, 0.0_dp), wscr(1:ngv, ngv0:ngv1,    imf, ipol), ngv  )
                 !
              END DO
              !
              CALL stop_clock("zgemm(make_scrn)")
              !
           END DO ! imf = 0, nmf
           !
        END DO ! ibnd = 1, nb(1)
        !
        ! Recieving WFCs from non-blocking circular shift
        !
        CALL mp_circular_shift_left(nb(2), 1, band2_comm )
        IF(okvan) THEN
           CALL mp_waitall(req(1:6))
        ELSE
           CALL mp_waitall(req(1:4))
        END IF
        wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 2) = wfcq_r(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz)
        IF(okvan) becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 2) &
        &     = becwfcq_r(1:nkb, 1:nb_max, 1:npol, 1:nqbz)
        wght(0:nmf, 1:nqbz, 1:nb_max, 1:nb_max) = wght_r(0:nmf, 1:nqbz, 1:nb_max, 1:nb_max)
        !
     END DO ! iproc2 = 1, nproc_band2
  END DO ! iproc1 = 1, nproc_band1
  !
  IF(nwscr == 2) &
  &  wscr(1:ngv, ngv0:ngv1, 0:nmf, 2) = wscr(1:ngv, ngv0:ngv1, 0:nmf, 1)
  !
  DEALLOCATE(rho1, rho2, wght, wfcq_r, wght_r, becwfcq_r)
  !
  CALL stop_clock("make_scrn")
  !
END SUBROUTINE make_scrn
!>
!> Calc. Kel
!>
SUBROUTINE make_Kel()
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum, mp_circular_shift_left, mp_barrier, mp_waitall
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE exx, ONLY : dfftt
  USE io_global, ONLY : stdout
  USE uspp, ONLY : nkb, okvan
  !
  USE sctk_val, ONLY : gq2, ngv, nmf, nb, nb_max, wfcq, becwfcq, &
  &                    ngv0, ngv1, Kel, Wscr, lsf, nqbz
  USE mp_kel, ONLY : nproc_band1, nproc_band2, band1_comm, band2_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, jbnd, imf, iproc1, iproc2, ipol, jpol, npol2, &
  &          nwscr, nKel, jbnd_ik, req(4)
  REAL(dp) :: RAM_Kel
  !
  COMPLEX(dp),ALLOCATABLE :: rho1(:,:,:), rho2(:,:), Kel0(:,:), Kel_temp(:,:,:,:), &
  &                          wfcq_r(:,:,:,:), becwfcq_r(:,:,:,:)
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
  RAM_Kel = REAL(nmf+2,dp)*REAL(nb_max**2,DP)*REAL(nqbz)*REAL(nKel,DP)
  CALL mp_sum(RAM_Kel, world_comm)
  WRITE(stdout,'(9x,"Total RAM for Kel : ",e10.2," GB")') RAM_Kel*8.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(Kel(0:nmf+1,nqbz,nb_max,nb_max,nKel))
  Kel(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 1:nKel) = 0.0_dp
  !
  ALLOCATE(rho1(ngv,nqbz*nb_max,npol2), rho2(ngv0:ngv1, nqbz*nb_max), &
  &        Kel0(nqbz*nb_max,nwscr), wfcq_r(dfftt%nnr, nb_max, npol, nqbz))
  IF(okvan) THEN
     ALLOCATE(becwfcq_r(nkb, nb_max, npol, nqbz))
  ELSE
     ALLOCATE(becwfcq_r(1, 1, 1, 1))
  END IF
  !
  DO iproc1 = 1, nproc_band1
     !
     CALL mp_circular_shift_left(nb(1), 1, band1_comm )
     CALL circular_shift_wrapper_c(dfftt%nnr, nb_max*npol*nqbz, band1_comm, &
     &                             wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 1))
     IF(okvan) CALL circular_shift_wrapper_c(nkb, nb_max*npol*nqbz, band1_comm, &
     &                                       becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 1))
     CALL circular_shift_wrapper_r(nmf+2, nqbz*nb_max*nb_max*nKel, band1_comm, &
     &                             Kel(0:nmf+1, 1:nqbz, 1:nb_max, 1:nb_max, 1:nKel))
     !
     DO iproc2 = 1, nproc_band2
        !
        ! Non-blocking circular shift : They are used in the next loop.
        !
        CALL circular_shift_wrapper_c_nb(dfftt%nnr, nb_max*npol*nqbz, band2_comm, &
        &                                wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 2), &
        &                              wfcq_r(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz), req(1:2))
        IF(okvan) CALL circular_shift_wrapper_c_nb(nkb, nb_max*npol*nqbz, band2_comm, &
        &                                          becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 2), &
        &                                        becwfcq_r(1:nkb, 1:nb_max, 1:npol, 1:nqbz), req(3:4))
        !
        DO ibnd = 1, nb(1)
           !
           CALL calc_rhog(npol2, ibnd, rho1)
           !
           ! Infinite frequency -> Bare Coulomb
           !
           !$OMP PARALLEL DEFAULT(NONE) &
           !$OMP & SHARED(nb,nqbz,ngv,rho1,gq2,omega,nmf,Kel,ibnd) &
           !$OMP & PRIVATE(jbnd_ik,jbnd,ik)
           !$OMP DO
           DO jbnd_ik = 1, nb(2)*nqbz
              !
              jbnd = 1 + (jbnd_ik-1) / nqbz
              ik   = 1 +  jbnd_ik-1 -  nqbz*(jbnd-1)
              !
              Kel(nmf + 1,ik,jbnd,ibnd,1) = SUM(REAL(rho1(1:ngv,jbnd_ik,1) &
              &                              * CONJG(rho1(1:ngv,jbnd_ik,1)), dp) &
              &                              /        gq2(1:ngv))  / omega
              !
           END DO ! jbnd_ik = 1, nb(2)*nqbz
           !$OMP END DO
           !$OMP END PARALLEL
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
                 CALL zgemm("T", "N", ngv1-ngv0+1, nqbz*nb(2), ngv, &
                 &          (1.0_dp, 0.0_dp), Wscr(1:ngv, ngv0:ngv1, imf,ipol), ngv, &
                 &                            rho1(1:ngv,            1:nqbz*nb_max, jpol), ngv, &
                 &          (0.0_dp, 0.0_dp), rho2(       ngv0:ngv1, 1:nqbz*nb_max), ngv1 - ngv0 + 1)
                 !
                 CALL stop_clock("zgemm(make_kel)")
                 !
                 !$OMP PARALLEL DEFAULT(NONE) &
                 !$OMP & SHARED(nb,nqbz,nb_max,ngv0,ngv1,rho1,rho2,Kel0,ipol,jpol) &
                 !$OMP & PRIVATE(jbnd_ik)
                 !$OMP DO
                 DO jbnd_ik = 1, nb(2)*nqbz
                    Kel0(jbnd_ik, ipol) = DOT_PRODUCT(rho1(ngv0:ngv1, jbnd_ik, jpol), &
                    &                                 rho2(ngv0:ngv1, jbnd_ik      )  )
                    !
                 END DO ! jbnd_ik = 1, nb(2)*nqbz
                 !$OMP END DO
                 !$OMP END PARALLEL
                 !
              END DO ! ipol = 1, nwscr
              !
              !$OMP PARALLEL DEFAULT(NONE) &
              !$OMP & SHARED(nb,nqbz,omega,nmf,Kel,ibnd,nwscr,Kel0,imf) &
              !$OMP & PRIVATE(jbnd_ik,jbnd,ik)
              !$OMP DO
              DO jbnd_ik = 1, nb(2)*nqbz
                 !
                 jbnd = 1 + (jbnd_ik-1) / nqbz
                 ik   = 1 +  jbnd_ik-1 -  nqbz*(jbnd-1)
                 !
                 Kel(imf,ik,   jbnd,ibnd,1) = Kel(imf,ik,jbnd,ibnd,1) + REAL(Kel0(jbnd_ik,1), dp) / omega
                 !
                 IF(nwscr == 2) THEN
                    Kel(imf,ik,jbnd,ibnd,2) = Kel(imf,ik,jbnd,ibnd,2) + REAL(Kel0(jbnd_ik,2), dp) / omega
                 ELSE IF(nwscr == 4) THEN
                    Kel(imf,ik,jbnd,ibnd,2) = Kel(imf,ik,jbnd,ibnd,2) + REAL(SUM(Kel0(jbnd_ik,2:4)), dp) / omega
                 END IF
                 !
              END DO ! jbnd_ik = 1, nb(2)*nqbz
              !$OMP END DO
              !$OMP END PARALLEL
              !
           END DO ! imf = 0, nmf
        END DO
        !
        ! Recieving WFCs from non-blocking circular shift
        !
        CALL mp_circular_shift_left(nb(2), 1, band2_comm )
        IF(okvan) THEN
           CALL mp_waitall(req(1:4))
        ELSE
           CALL mp_waitall(req(1:2))
        END IF
        wfcq(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz, 2) = wfcq_r(1:dfftt%nnr, 1:nb_max, 1:npol, 1:nqbz)
        IF(okvan) becwfcq(1:nkb, 1:nb_max, 1:npol, 1:nqbz, 2) &
        &     = becwfcq_r(1:nkb, 1:nb_max, 1:npol, 1:nqbz)
        CALL circular_shift_wrapper_r(nmf+2, nqbz*nb_max*nb_max*nKel, band2_comm, &
        &                             Kel(0:nmf+1, 1:nqbz, 1:nb_max, 1:nb_max, 1:nKel))
        !
     END DO
  END DO
  !
  DEALLOCATE(rho1, rho2, Kel0, wfcq_r, becwfcq_r)
  !
  ! Average the weights for degenerated states
  !
  ALLOCATE(Kel_temp(0:nmf + 1,nqbz,nb_max,nb_max))
  IF(lsf>0) THEN
     Kel_temp(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max) = CMPLX(Kel(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 1), &
     &                                                       Kel(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 2), DP)
  ELSE
     Kel_temp(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max) = CMPLX(Kel(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 1), 0.0_dp, DP)
  END IF
  CALL average_wght_kel(nmf+2, Kel_temp(0:nmf+1,1:nqbz,1:nb_max,1:nb_max))
  Kel(          0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 1) = REAL( Kel_temp(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max), DP)
  IF(lsf>0) Kel(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max, 2) = AIMAG(Kel_temp(0:nmf + 1, 1:nqbz, 1:nb_max, 1:nb_max))
  !
  DEALLOCATE(Kel_temp)
  !
  CALL stop_clock("make_kel")
  !
END SUBROUTINE make_Kel
!>
!>
!>
SUBROUTINE average_wght_kel(nmf0, wght)
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  !
  USE sctk_val, ONLY : degen, ndegen, nqbz, nb, bdsp, nb_max, kplusq
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nmf0
  COMPLEX(DP),INTENT(INOUT) :: wght(nmf0,nqbz, nb_max, nb_max)
  !
  INTEGER :: ik, jk, ibnd, idegen
  COMPLEX(DP) :: wght_ave(nmf0, nbnd)
  COMPLEX(DP),ALLOCATABLE :: wght_tmp(:,:,:)
  !
  ALLOCATE(wght_tmp(nmf0,nbnd,nbnd))
  !
  DO ik = 1, nqbz
     !
     wght_tmp(1:nmf0, 1:nbnd, 1:nbnd) = 0.0_dp
     wght_tmp(1:nmf0,    bdsp(2)+1:bdsp(2)+nb(2), bdsp(1)+1:bdsp(1)+nb(1)) &
     & = wght(1:nmf0, ik,        1:        nb(2),         1:        nb(1))
     CALL mp_sum(wght_tmp, world_comm)
     !
     DO idegen = 1, ndegen(ik, 1)
        wght_ave(1:nmf0, 1:nbnd) = SUM(wght_tmp(1:nmf0, 1:nbnd, degen(1, idegen, ik, 1): &
        &                                                       degen(2, idegen, ik, 1)), 3) &
        &                      / REAL(degen(2, idegen, ik,1) - degen(1,idegen,ik,1) + 1, DP)
        !
        DO ibnd = degen(1,idegen, ik, 1), degen(2, idegen, ik, 1)
           wght_tmp(1:nmf0, 1:nbnd, ibnd)  = wght_ave(1:nmf0, 1:nbnd)
        END DO
     END DO
     !
     jk = kplusq(ik)
     DO idegen = 1, ndegen(jk, 2)
        wght_ave(1:nmf0,1 :nbnd) = SUM(wght_tmp(1:nmf0, degen(1, idegen, jk, 2): &
        &                                               degen(2, idegen, jk, 2), 1:nbnd), 2) &
        &                      / REAL(degen(2, idegen, jk, 2) - degen(1, idegen, jk, 2) + 1, DP)
        !
        DO ibnd = degen(1, idegen, jk, 2), degen(2, idegen, jk, 2)
           wght_tmp(1:nmf0, ibnd, 1:nbnd)  = wght_ave(1:nmf0, 1:nbnd)
        END DO
     END DO
     !
     wght(        1:nmf0, ik,     1:        nb(2),         1:        nb(1)) &
     & = wght_tmp(1:nmf0, bdsp(2)+1:bdsp(2)+nb(2), bdsp(1)+1:bdsp(1)+nb(1))
     !
  END DO
  !
  DEALLOCATE(wght_tmp)
  !
END SUBROUTINE average_wght_kel
!>
!>
!>
SUBROUTINE calc_rhog(npol2, ibnd, rho1)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE fft_scalar, ONLY : cfft3d
  USE noncollin_module, ONLY : npol
  USE exx, ONLY : dfftt
  USE uspp, ONLY : nkb, okvan
  !
  USE sctk_val, ONLY : gindx, ngv, wfcq, becwfcq, nqbz, nb_max, nb
  !
  IMPLICIT none
  !
  INTEGER,INTENT(IN) :: npol2, ibnd
  COMPLEX(DP),INTENT(OUT) :: rho1(ngv, nqbz*nb_max, npol2)
  !
  INTEGER:: jbnd, ik, ipol, jbnd_ik
  COMPLEX(dp),ALLOCATABLE :: rho0(:,:,:), rho4(:)
  !
  ALLOCATE(rho0(dfftt%nnr,nqbz*nb_max, npol2))
  !
  CALL start_clock("wfc*wfc=rho")
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nb,nqbz,npol2,dfftt,npol,wfcq,becwfcq,omega,ibnd,okvan,rho0,nkb) &
  !$OMP & PRIVATE(jbnd_ik, jbnd, ik, ipol, rho4)
  !
  ALLOCATE(rho4(dfftt%nnr))
  !
  !$OMP DO
  DO jbnd_ik = 1, nb(2)*nqbz
     !
     jbnd = 1 + (jbnd_ik-1) / nqbz
     ik   = 1 +  jbnd_ik-1 -  nqbz*(jbnd-1)
     !
     ! debug
     !rho4(1:dfftt%nnr) = conjg(wfcq(1:dfftt%nnr,ibnd,1,ik, 1)) * wfcq(1:dfftt%nnr,ibnd,1,ik, 1)
     !write(*,'(a,2f15.8)',advance="no") "debug1", sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
     !IF(okvan) then
     !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) / omega
     !   CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd,1,ik,1), &
     !   &                                 becwfcq(1:nkb,ibnd,1,ik,1))
     !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) * omega
     !   write(*,'(2f15.8)',advance="no") sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
     !END IF
     !rho4(1:dfftt%nnr) = conjg(wfcq(1:dfftt%nnr,ibnd,1,ik, 2)) * wfcq(1:dfftt%nnr,ibnd,1,ik, 2)
     !write(*,'(a,2f15.8)',advance="no") "debug1", sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
     !IF(okvan) then
     !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) / omega
     !   CALL addusxx_r(rho4(1:dfftt%nnr), becwfcq(1:nkb,ibnd,1,ik,2), &
     !   &                                 becwfcq(1:nkb,ibnd,1,ik,2))
     !   rho4(1:dfftt%nnr) = rho4(1:dfftt%nnr) * omega
     !   write(*,'(2f15.8)') sum(rho4(1:dfftt%nnr)) / dble(dfftt%nnr)
     !END IF
     ! end debug
     !
     rho0(1:dfftt%nnr, jbnd_ik, 1:npol2) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
     DO ipol = 1, npol
        !
        rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr, ibnd, ipol, ik, 1)) &
        &                       * wfcq(1:dfftt%nnr, jbnd, ipol, ik, 2) / omega
        IF(okvan) CALL addus_kel(rho4(1:dfftt%nnr), becwfcq(1:nkb, ibnd, ipol, ik, 1), &
        &                                           becwfcq(1:nkb, jbnd, ipol, ik, 2))
        rho0(1:dfftt%nnr, jbnd_ik,1) = rho0(1:dfftt%nnr, jbnd_ik,1) + rho4(1:dfftt%nnr) * omega
        !
        IF(npol2 == 4) THEN
           !
           rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr, ibnd,   ipol, ik, 1)) &
           &                       * wfcq(1:dfftt%nnr, jbnd, 3-ipol, ik, 2) / omega
           IF(okvan) CALL addus_kel(rho4(1:dfftt%nnr), becwfcq(1:nkb, ibnd,   ipol, ik, 1), &
           &                                           becwfcq(1:nkb, jbnd, 3-ipol, ik, 2))
           rho0(1:dfftt%nnr, jbnd_ik,2) = rho0(1:dfftt%nnr, jbnd_ik,2) + rho4(1:dfftt%nnr) * omega
           !
           rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr, ibnd,   ipol, ik, 1)) &
           &                       * wfcq(1:dfftt%nnr, jbnd, 3-ipol, ik, 2) / omega
           IF(okvan) CALL addus_kel(rho4(1:dfftt%nnr), becwfcq(1:nkb, ibnd,   ipol, ik, 1), &
           &                                           becwfcq(1:nkb, jbnd, 3-ipol, ik, 2))
           rho0(1:dfftt%nnr, jbnd_ik,3) = rho0(1:dfftt%nnr, jbnd_ik,3) &
           &                   + rho4(1:dfftt%nnr) * REAL(3-2*ipol, dp) * omega
           !
           rho4(1:dfftt%nnr) = CONJG(wfcq(1:dfftt%nnr, ibnd, ipol, ik, 1)) &
           &                       * wfcq(1:dfftt%nnr, jbnd, ipol, ik, 2) / omega
           IF(okvan) CALL addus_kel(rho4(1:dfftt%nnr), becwfcq(1:nkb, ibnd, ipol, ik, 1), &
           &                                           becwfcq(1:nkb, jbnd, ipol, ik, 2))
           rho0(1:dfftt%nnr, jbnd_ik,4) = rho0(1:dfftt%nnr, jbnd_ik,4) &
           &                   + rho4(1:dfftt%nnr) * REAL(3-2*ipol, dp) * omega
           !
        END IF
        !
     END DO ! ipol
     !
  END DO ! jbnd_ik = nb(2)*nqbz
  !$OMP END DO
  !
  DEALLOCATE(rho4)
  !
  !$OMP END PARALLEL
  !
  CALL stop_clock("wfc*wfc=rho")
  !
  CALL start_clock("fft[rho]")
  !
  DO ipol = 1, npol2
     DO jbnd_ik = 1, nb(2)*nqbz
        !
        CALL cfft3d (rho0(1:dfftt%nnr, jbnd_ik, ipol), &
        & dfftt%nr1, dfftt%nr2, dfftt%nr3, dfftt%nr1, dfftt%nr2, dfftt%nr3, 1, -1)
        !
        rho1(1:ngv, jbnd_ik, ipol) = rho0(gindx(1:ngv), jbnd_ik, ipol)
        !
     END DO ! jbnd_ik = nb(2)*nqbz
  END DO ! ipol
  !
  CALL stop_clock("fft[rho]")
  !
  DEALLOCATE(rho0)
  !
END SUBROUTINE calc_rhog
  !
  !
  !------------------------------------------------------------------------
SUBROUTINE addus_kel( rho, becphi, becpsi )
   !------------------------------------------------------------------------
   !! This routine adds to the two wavefunctions density (in real space) 
   !! the part that is due to the US augmentation.  
   !! NOTE: the density in this case is NOT real and NOT normalized to 1, 
   !!       except when (bec-)\(\text{phi}\) and (bec-)\(\text{psi}\) are equal,
   !!       or with gamma tricks.
   !
   !! With gamma tricks: input rho must contain contributions from band 1
   !! in real part, from band 2 in imaginary part. Call routine twice:
   !
   !! * with \(\text{becphi}=\langle\beta|\phi(1)\rangle\) (real);
   !! * then with \(\text{becphi}=-i\langle\beta|\phi(2)\rangle\) (imaginary).
   !
   USE kinds,            ONLY : DP
   USE ions_base,        ONLY : nat, ityp
   USE uspp,             ONLY : nkb, ijtoh, ofsbeta
   USE uspp_param,       ONLY : upf, nh
   USE realus,           ONLY : tabxx
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(INOUT) :: rho(:)
   !! charge density
   COMPLEX(DP), INTENT(IN) :: becphi(nkb)
   !! see main comment
   COMPLEX(DP), INTENT(IN) :: becpsi(nkb)
   !! see main comment
   !
   ! ... local variables
   !
   INTEGER :: ia, nt, ir, irb, ih, jh, mbia
   INTEGER :: ikb, jkb
   !
   DO ia = 1, nat
     !
     mbia = tabxx(ia)%maxbox
     IF ( mbia == 0 ) CYCLE
     !
     nt = ityp(ia)
     IF ( .NOT. upf(nt)%tvanp ) CYCLE
     !
     DO ih = 1, nh(nt)
       DO jh = 1, nh(nt)
         ikb = ofsbeta(ia) + ih
         jkb = ofsbeta(ia) + jh
         DO ir = 1, mbia
           irb = tabxx(ia)%box(ir)
           rho(irb) = rho(irb) + tabxx(ia)%qr(ir,ijtoh(ih,jh,nt)) &
                                 *CONJG(becphi(ikb))*becpsi(jkb)
         ENDDO
       ENDDO
     ENDDO
     !
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE addus_kel
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
  USE sctk_val, ONLY : nmf, mf, Kel, lsf, nb, bdsp
  !
  IMPLICIT NONE
  !
  INTEGER :: ntetra0, ntetra1, ii, ikv(3), ik, imf, kintp(8), nkel, ikel, &
  &          b_low_p(2), b_high_p(2), b_low_g(2), b_high_g(2)
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
  DO ii = 1, 2
     b_low_g( ii) = MAX(b_low,  bdsp(ii)+1)
     b_high_g(ii) = MIN(b_high, bdsp(ii)+nb(ii))
     b_low_p( ii) = b_low_g( ii) - bdsp(ii)
     b_high_p(ii) = b_high_g(ii) - bdsp(ii)
  END DO
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks,nk1,nk2,nk3,nq1,nq2,nq3,nkel,nmf,Kel,b_low_p,b_high_p,b_low_g,b_high_g,wght,mu) &
  !$OMP & PRIVATE(ik,ikv,kv,ikel,imf,ii,kintp,wintp)
  !$OMP DO REDUCTION(+:mu)
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
              mu(imf,ikel) = mu(imf,ikel) + &
              & SUM(Kel(imf, kintp(ii), b_low_p(2):b_high_p(2), b_low_p(1):b_high_p(1), ikel) &
              &     * wght(             b_low_g(2):b_high_g(2), b_low_g(1):b_high_g(1), ik)) * wintp(ii)
              !
           END DO ! ii = 1, 8
        END DO ! imf = 0, nmf + 1
     END DO ! ikel = 1, nkel
     !
  END DO ! ik = 1, nks
  !$OMP END DO
  !$OMP END PARALLEL
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
  USE sctk_val, ONLY : Kel, nqbz, nmf, lsf, mf, nci, nqbz, nb
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jb, imf, jmf, ikel, nkel, ik_jb_ib
  REAL(dp) :: coef(nci,nci), Kel0(nci), Kel1(0:nmf+1), cheb(1:nci,0:nmf+1), &
  &           err0(0:nmf+1), aveer(0:nmf+1), maxer(0:nmf+1), mf0, x0
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
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(nqbz,nb,nci,ikel,coef,Kel,cheb,nmf,aveer,maxer) &
     !$OMP & PRIVATE(ik_jb_ib,ib,jb,ik,Kel1,Kel0,imf,jmf,err0)
     !$OMP DO REDUCTION(+:aveer) REDUCTION(MAX:maxer)
     DO ik_jb_ib = 1, nqbz*nb(2)*nb(1)
        !
        ib = 1 + (ik_jb_ib-1) / (nqbz*nb(2))
        jb = 1 + (ik_jb_ib-1  -  nqbz*nb(2)*(ib-1)) / nqbz
        ik = 1 +  ik_jb_ib-1  -  nqbz*nb(2)*(ib-1)  - nqbz*(jb-1)
        !
        Kel0(1:nci) = 0.0_dp
        DO imf = 1, nci
           DO jmf = 1, nci
              Kel0(imf) = Kel0(imf) + coef(imf,jmf)*Kel((jmf-1)*2,ik,jb,ib,ikel)
           END DO
        END DO
        !
        Kel1(0:nmf+1) = MATMUL(Kel0(1:nci), cheb(1:nci, 0:nmf+1))
        !
        err0(0:nmf+1) = ABS((Kel1(0:nmf+1)  -     Kel(0:nmf+1,ik,jb,ib,ikel))) &
        &   * 2.0_dp / (ABS( Kel1(0:nmf+1)) + ABS(Kel(0:nmf+1,ik,jb,ib,ikel)))
        !
        aveer(0:nmf+1) = aveer(0:nmf+1) + err0(0:nmf+1)
        maxer(0:nmf+1) = MAX(maxer(0:nmf+1), err0(0:nmf+1))
        !
        Kel(0:nci-1,ik,jb,ib,ikel) = Kel0(1:nci)
        !
     END DO ! ik_jb_ib = 1, nqbz*nb(1)*nb(2)
     !$OMP END DO NOWAIT
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
  USE sctk_val, ONLY : Kel, nqbz, lsf, nci, nb, bdsp, nmf, nb_max
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, nkel, iproc, ik
  CHARACTER(6), EXTERNAL :: int_to_char
  REAL(DP),ALLOCATABLE :: Kel_all(:,:,:,:,:)
  !
  IF(lsf>0) THEN
     nkel = 2
  ELSE
     nkel = 1
  END IF
  IF(mpime==0) THEN
     ALLOCATE(Kel_all(0:nci-1,nbnd,nbnd,nqbz,nkel))
  ELSE
     ALLOCATE(Kel_all(1,1,1,1,1))
  END IF
  !
  DO iproc = 1, nproc
     CALL mp_circular_shift_left(nb, 1, world_comm )
     CALL mp_circular_shift_left(bdsp, 1, world_comm )
     CALL circular_shift_wrapper_r(nmf + 2,nqbz*nb_max*nb_max*nkel, world_comm, Kel)
     DO ik = 1, nqbz
        IF(mpime==0) Kel_all(0:nci-1,     bdsp(2)+1:bdsp(2)+nb(2), bdsp(1)+1:bdsp(1)+nb(1), ik,1:nkel) &
        &              = Kel(0:nci-1, ik,         1:        nb(2),         1:        nb(1),    1:nkel)
        !
     END DO
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
  DEALLOCATE(Kel_all, Kel)
  !
END SUBROUTINE write_Kel
!
END MODULE sctk_coulomb
