!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_wfc
  !
  IMPLICIT NONE
  !
CONTAINS
!
! read g-vector & wfc
!
SUBROUTINE get_wfcg()
  !
  USE wvfct, ONLY : nbnd, npwx
  USE kinds, ONLY : DP
  USE io_files, ONLY : prefix, tmp_dir
  USE mp_world, ONLY : world_comm
  USE mp_pools, ONLY : intra_pool_comm, inter_pool_comm, me_pool, nproc_pool
  USE mp, ONLY : mp_max, mp_sum, mp_bcast, mp_circular_shift_left, mp_barrier
  USE disp, ONLY : nq1, nq2, nq3
  USE io_global, ONLY : stdout
  USE noncollin_module, ONLY : npol
  !
  USE sctk_coulomb, ONLY : circular_shift_wrapper_c
  USE sctk_val, ONLY : nqbz, wfc0, igv, npw, nk_p, nk_p_max, k0_p, nbnd_p, nbnd_p_max, bnd0_p
  !
  IMPLICIT NONE
  !
  INTEGER :: fi, ii, ngw0
  CHARACTER(LEN=320) :: filename
  !
  INTEGER :: ik0, ispin, npol0, ierr, ikv(3), ibnd, iproc, ik_p, ik_g
  REAL(dp) :: scalef, xk0(3), bvec0(3,3), xk1(3), RAM_wfc
  LOGICAL :: gamma_only
  COMPLEX(DP),ALLOCATABLE :: wfc0_all(:,:,:,:,:)
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  fi = find_free_unit()
  !
  CALL divide(inter_pool_comm, nqbz, k0_p, nk_p)
  k0_p = k0_p - 1
  nk_p = nk_p - k0_p
  nk_p_max = nk_p
  CALL mp_max(nk_p_max, world_comm)
  WRITE(stdout,'(7x,"k-point per PE : ",i0)') nk_p_max
  !
  CALL divide(intra_pool_comm, nbnd, bnd0_p, nbnd_p)
  bnd0_p = bnd0_p - 1
  nbnd_p = nbnd_p - bnd0_p
  nbnd_p_max = nbnd_p
  CALL mp_max(nbnd_p_max, world_comm)
  WRITE(stdout,'(7x,"Bands per PE : ",i0)') nbnd_p_max
  !
  CALL mp_sum(npwx, intra_pool_comm)
  CALL mp_max(npwx, inter_pool_comm)
  !
  RAM_wfc = REAL(npwx,dp)*REAL(nbnd_p_max,dp)*REAL(npol,dp)*REAL(nk_p_max,dp)*2.0_dp
  CALL mp_sum(RAM_wfc, world_comm)
  WRITE(stdout,'(7x,"Total RAM for input WFC : ",e10.2," GB")') RAM_wfc*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wfc0(npwx,nbnd_p_max,npol,nk_p_max,2), igv(3,npwx,nk_p_max,2), &
  &        npw(nk_p_max,2))
  igv(1:3,1:npwx,                    1:nk_p_max,1:2) = 0
  wfc0(   1:npwx,1:nbnd_p_max,1:npol,1:nk_p_max,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Read Wfc(k)
  !
  IF(me_pool == 0) THEN
     !
     ALLOCATE(wfc0_all(npwx,nbnd,npol,nk_p_max,2))
     !
     DO ii = 1, 2
        !
        DO ik_p = 1, nk_p
           !
           ik_g = ik_p + k0_p
           !
           filename = TRIM(tmp_dir) // TRIM(prefix) // '.save/' // 'wfc' // &
           &          TRIM(int_to_char(nqbz * (ii-1) + ik_g))
           !
           OPEN(unit=fi, file=trim(filename)//'.dat', form='unformatted', status='old', iostat=ierr)
           READ(fi) ik0, xk0(1:3), ispin, gamma_only, scalef
           READ(fi) ngw0, npw(ik_p,ii), npol0, ibnd
           READ(fi) bvec0(1:3,1:3)
           READ(fi) igv(1:3,1:npw(ik_p,ii),ik_p,ii)
           DO ibnd = 1, nbnd
              READ(fi) wfc0_all(1:npw(ik_p,ii),ibnd,1:npol,ik_p,ii)
           END DO
           CLOSE(unit=fi, status='keep')
           !
           ! Verify k-point
           !
           ikv(1) = (ik_g - 1) / (nq3*nq2)
           ikv(2) = (ik_g - 1 - ikv(1)*nq2*nq3) / nq3
           ikv(3) =  ik_g - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
           WHERE(ikv(1:3)*2 + ii - 1 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
           !
           xk1(1:3) = (REAL(ikv(1:3), dp) + 0.5_dp*REAL(ii-1, dp)) / REAL((/nq1, nq2, nq3/), dp)
           xk1(1:3) = MATMUL(bvec0(1:3, 1:3), xk1(1:3))
           IF(ANY(ABS(xk1(1:3) - xk0(1:3)) > 1.0e-5_dp)) THEN
              WRITE(stdout,'("ik = ",i5,", ii = ",i5)') ik_g, ii
              WRITE(stdout,'("In ",3f15.10,", Required ",3f15.10)') xk0(1:3), xk1(1:3)
              CALL errore ('get_wfcg', 'Invalid k-point (order)', ik_g)
           END IF
           !
        END DO ! ik_p
        !
     END DO ! ii
     !
  END IF
  !
  ! Scatter wfc0_all(bnd0_p+1:bnd0_p+nbnd_p) -> wfc0(1:nbnd_p) intra pool
  !
  CALL mp_bcast(npw,0,intra_pool_comm)
  CALL mp_bcast(igv,0,intra_pool_comm)
  !
  DO iproc = 1, nproc_pool
     IF(me_pool == 0) wfc0(1:npwx,       1:       nbnd_p,1:npol,1:nk_p_max,1:2) &
     &          = wfc0_all(1:npwx,bnd0_p+1:bnd0_p+nbnd_p,1:npol,1:nk_p_max,1:2)
     CALL circular_shift_wrapper_c(npwx,nbnd_p_max*npol*nk_p_max*2,intra_pool_comm,wfc0)
     CALL mp_circular_shift_left(nbnd_p, 1, intra_pool_comm)
     CALL mp_circular_shift_left(bnd0_p, 1, intra_pool_comm)
  END DO
  !
  IF(me_pool == 0) DEALLOCATE(wfc0_all)
  !
END SUBROUTINE get_wfcg
!
! FFT wavefunctions
!
SUBROUTINE fft_wfc()
  !
  USE wvfct, ONLY : nbnd, npwx
  USE kinds, ONLY : DP
  USE mp, ONLY : mp_max, mp_sum
  USE mp_world, ONLY : world_comm
  USE fft_scalar, ONLY : cfft3d
  USE io_global, ONLY : stdout
  USE noncollin_module, ONLY : npol
  USE exx, ONLY : dfftt
  USE gvect, ONLY : ngm, mill
  USE becmod, ONLY : calbec
  USE cell_base, ONLY : bg
  USE uspp, ONLY : nkb, okvan
  USE disp, ONLY : nq1, nq2, nq3
  !
  USE sctk_val, ONLY : igv, npw, wfc0, k0_p, nk_p, nk_p_max, nbnd_p, nbnd_p_max, &
  &                    wfc, wfcq, becwfc, becwfcq, nqbz
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ik_p, ik_g, ib, igv2(3), ig2, ipol, ikv(3), ii, igm
  REAL(DP) :: xk(3), RAM_wfc
  INTEGER,ALLOCATABLE :: igk(:)
  COMPLEX(DP),ALLOCATABLE :: vkb(:,:)
  !
  ALLOCATE(vkb(npwx,nkb), igk(npwx))
  !
  RAM_wfc = REAL(dfftt%nnr, dp)*REAL(nbnd_p_max, dp)*REAL(npol, dp)*REAL(nk_p_max, dp)*4.0_dp
  CALL mp_sum(RAM_wfc, world_comm)
  WRITE(stdout,'(7x,"Total RAM for WFC(r) : ",e10.2," GB")') &
  &     REAL(dfftt%nnr,dp)*REAL(npol*4,dp)*REAL(nbnd,dp)*REAL(nqbz,dp)*16.0e-9_dp
  !
  ALLOCATE(wfc( dfftt%nnr, nbnd_p_max, npol, nk_p_max, 2), &
  &        wfcq(dfftt%nnr, nbnd_p_max, npol, nk_p_max, 2))
  ALLOCATE(becwfc(    nkb, nbnd_p_max, npol, nk_p_max, 2), &
  &        becwfcq(   nkb, nbnd_p_max, npol, nk_p_max, 2))
  wfc(       1:dfftt%nnr,1:nbnd_p_max,1:npol,1:nk_p_max,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  IF(okvan) becwfc(1:nkb,1:nbnd_p_max,1:npol,1:nk_p_max,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Fourier trans. wfc(G) -> wfc(r)
  !
  DO ii = 1, 2
     DO ik_p = 1, nk_p
        !
        DO ig = 1, npw(ik_p,ii)
           igv2(1:3) = MODULO(igv(1:3,ig,ik_p,ii), (/dfftt%nr1, dfftt%nr2, dfftt%nr3/))
           ig2 = 1 + igv2(1) + dfftt%nr1 * igv2(2) + dfftt%nr1 * dfftt%nr2 * igv2(3)
           wfc(ig2, 1:nbnd_p, 1:npol, ik_p, ii) = wfc0(ig, 1:nbnd_p, 1:npol, ik_p, ii)
        END DO
        !
        DO ipol = 1, npol
           DO ib = 1, nbnd_p
              CALL cfft3d (wfc(1:dfftt%nnr, ib, ipol, ik_p, ii), &
              & dfftt%nr1, dfftt%nr2, dfftt%nr3, dfftt%nr1, dfftt%nr2, dfftt%nr3, 1, 1)
           END DO
        END DO
        !
     END DO
  END DO
  !
  ! <beta|phi_k>
  !
  IF ( okvan ) THEN
     DO ii = 1, 2
        DO ik_p = 1, nk_p
           !
           ik_g = ik_p + k0_p
           !
           ikv(1) = (ik_g - 1) / (nq3*nq2)
           ikv(2) = (ik_g - 1 - ikv(1)*nq2*nq3) / nq3
           ikv(3) =  ik_g - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
           WHERE(ikv(1:3)*2 + ii - 1 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
           !
           xk(1:3) = (REAL(ikv(1:3), dp) + 0.5_dp*REAL(ii-1, dp)) / REAL((/nq1, nq2, nq3/), dp)
           xk(1:3) = MATMUL(bg(1:3, 1:3), xk(1:3))
           !
           DO ig = 1, npw(ik_p,ii)
              DO igm = 1, ngm
                 IF(ALL(igv(1:3,ig,ik_p,ii) == mill(1:3,igm))) THEN
                   igk(ig) = igm
                   EXIT
                 END IF
              END DO
           END DO
           CALL init_us_2(npw(ik_p,ii), igk, xk, vkb)
           DO ipol = 1, npol
              CALL calbec(npw(ik_p,ii), vkb, wfc0(1:npwx, 1:nbnd_p, ipol,ik_p,ii), &
              &           becwfc(1:nkb, 1:nbnd_p, ipol, ik_p, ii), nbnd_p)
           END DO
        END DO
     END DO
  END IF
  !
  DEALLOCATE(wfc0, igv, npw, igk, vkb)
  !
END SUBROUTINE fft_wfc
!
END MODULE sctk_wfc
