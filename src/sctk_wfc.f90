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
  USE mp_world, ONLY : world_comm, mpime
  USE mp_pools, ONLY : inter_pool_comm, intra_pool_comm
  USE mp, ONLY : mp_max, mp_sum, mp_bcast, mp_circular_shift_left, mp_barrier
  USE disp, ONLY : nq1, nq2, nq3
  USE io_global, ONLY : stdout
  USE noncollin_module, ONLY : npol
  !
  USE sctk_coulomb, ONLY : circular_shift_wrapper_c
  USE sctk_val, ONLY : nqbz, wfc0, igv, npw, nb, bdsp, nb_max
  USE mp_kel, ONLY : mp_start_kel, band1_comm, band2_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: fi, ii, ngw0
  CHARACTER(320) :: filename
  !
  INTEGER :: ik0, ispin, npol0, ierr, ikv(3), ibnd, ik, npwx0
  REAL(dp) :: scalef, xk0(3), bvec0(3,3), xk1(3), RAM_wfc
  LOGICAL :: gamma_only
  COMPLEX(DP),ALLOCATABLE :: wfc0k(:,:,:)
  !
  CHARACTER(6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  fi = find_free_unit()
  !
  CALL mp_start_kel()
  CALL divide(band1_comm, nbnd, bdsp(1), nb(1))
  CALL divide(band2_comm, nbnd, bdsp(2), nb(2))
  !
  bdsp(1) = bdsp(1) - 1
  nb(1) = nb(1) - bdsp(1)
  !
  bdsp(2) = bdsp(2) - 1
  nb(2) = nb(2) - bdsp(2)
  !
  nb_max = MAXVAL(nb(1:2))
  CALL mp_max(nb_max, world_comm)
  !
  ! Compute maximum number of PWs across all k-points
  !
  npwx = 0
  IF(mpime == 0) THEN
    DO ik = 1, nqbz*2
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.save/' // 'wfc' // TRIM(int_to_char(ik))
      OPEN(unit=fi, file=trim(filename)//'.dat', form='unformatted', status='old', iostat=ierr)
      READ(fi) ik0, xk0(1:3), ispin, gamma_only, scalef
      READ(fi) ngw0, npwx0, npol0, ibnd
      CLOSE(unit=fi, status='keep')
      !
      npwx = MAX(npwx, npwx0)
    END DO
  END IF
  CALL mp_bcast(npwx, 0, world_comm)
  !
  RAM_wfc = REAL(npwx,dp)*REAL(nb_max,dp)*REAL(npol,dp)*REAL(nqbz,dp)*2.0_dp
  CALL mp_sum(RAM_wfc, world_comm)
  WRITE(stdout,'(7x,"Total RAM for input WFC : ",e10.2," GB")') RAM_wfc*16.0e-9_dp
  CALL mp_barrier(world_comm)
  ALLOCATE(wfc0(npwx,nb_max,npol,nqbz,2), igv(3,npwx,nqbz,2), &
  &        npw(nqbz,2), wfc0k(npwx,nbnd,npol))
  igv(1:3,1:npwx,                    1:nqbz,1:2) = 0
  wfc0(   1:npwx,1:nb_max,1:npol,1:nqbz,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Read Wfc(k)
  !
  DO ii = 1, 2
     !
     DO ik = 1, nqbz
        !
        IF(mpime == 0) THEN
           !
           filename = TRIM(tmp_dir) // TRIM(prefix) // '.save/' // 'wfc' // &
           &          TRIM(int_to_char(nqbz * (ii-1) + ik))
           !
           OPEN(unit=fi, file=trim(filename)//'.dat', form='unformatted', status='old', iostat=ierr)
           READ(fi) ik0, xk0(1:3), ispin, gamma_only, scalef
           READ(fi) ngw0, npw(ik,ii), npol0, ibnd
           READ(fi) bvec0(1:3,1:3)
           READ(fi) igv(1:3,1:npw(ik,ii),ik,ii)
           DO ibnd = 1, nbnd
              READ(fi) wfc0k(1:npw(ik,ii),ibnd,1:npol)
           END DO
           CLOSE(unit=fi, status='keep')
           !
           ! Verify k-point
           !
           ikv(1) = (ik - 1) / (nq3*nq2)
           ikv(2) = (ik - 1 - ikv(1)*nq2*nq3) / nq3
           ikv(3) =  ik - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
           WHERE(ikv(1:3)*2 + ii - 1 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
           !
           xk1(1:3) = (REAL(ikv(1:3), dp) + 0.5_dp*REAL(ii-1, dp)) / REAL((/nq1, nq2, nq3/), dp)
           xk1(1:3) = MATMUL(bvec0(1:3, 1:3), xk1(1:3))
           IF(ANY(ABS(xk1(1:3) - xk0(1:3)) > 1.0e-5_dp)) THEN
              WRITE(stdout,'("ik = ",i5,", ii = ",i5)') ik, ii
              WRITE(stdout,'("In ",3f15.10,", Required ",3f15.10)') xk0(1:3), xk1(1:3)
              CALL errore ('get_wfcg', 'Invalid k-point (order)', ik)
           END IF
           !
        END IF
        !
        CALL mp_bcast(npw(            ik, ii),         0, world_comm)
        CALL mp_bcast(igv(1:3,1:npwx, ik, ii),         0, world_comm)
        CALL mp_bcast(wfc0k  (1:npwx, 1:nbnd, 1:npol), 0, world_comm)
        !
        wfc0(1:npwx, 1:nb(ii),1:npol,ik,ii) = wfc0k(1:npwx,bdsp(ii)+1:bdsp(ii)+nb(ii),1:npol)
        !
     END DO ! ik
     !
  END DO ! ii
  !
  DEALLOCATE(wfc0k)
  !
END SUBROUTINE get_wfcg
!
! FFT wavefunctions
!
SUBROUTINE fft_wfc()
  !
  USE wvfct, ONLY : npwx
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
  USE uspp_init, ONLY : init_us_2
  !
  USE sctk_val, ONLY : igv, npw, wfc0, nb, nb_max, &
  &                    wfc, wfcq, becwfc, becwfcq, nqbz
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ik, ib, igv2(3), ig2, ipol, ikv(3), ii, igm, ii_ik
  REAL(DP) :: xk(3), RAM_wfc
  INTEGER,ALLOCATABLE :: igk(:)
  COMPLEX(DP),ALLOCATABLE :: vkb(:,:)
  !
  ALLOCATE(vkb(npwx,nkb), igk(npwx))
  !
  RAM_wfc = REAL(dfftt%nnr, dp)*REAL(nb_max, dp)*REAL(npol, dp)*REAL(nqbz, dp)*4.0_dp
  CALL mp_sum(RAM_wfc, world_comm)
  WRITE(stdout,'(7x,"Total RAM for WFC(r) : ",e10.2," GB")') RAM_wfc*16.0e-9_dp
  !
  ALLOCATE(wfc( dfftt%nnr, nb_max, npol, nqbz, 2), &
  &        wfcq(dfftt%nnr, nb_max, npol, nqbz, 2))
  wfc(       1:dfftt%nnr,1:nb_max,1:npol,1:nqbz,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  IF(okvan) THEN
     !
     ALLOCATE(becwfc(    nkb, nb_max, npol, nqbz, 2), &
     &        becwfcq(   nkb, nb_max, npol, nqbz, 2))
     becwfc(1:nkb,1:nb_max,1:npol,1:nqbz,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
     !
     RAM_wfc = REAL(nkb, dp)*REAL(nb_max, dp)*REAL(npol, dp)*REAL(nqbz, dp)*4.0_dp
     CALL mp_sum(RAM_wfc, world_comm)
     WRITE(stdout,'(7x,"Total RAM for <beta|WFC> : ",e10.2," GB")') RAM_wfc*16.0e-9_dp
     !
  END IF
  !
  ! Fourier trans. wfc(G) -> wfc(r)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nqbz,npw,dfftt,nb,npol,wfc,wfc0,igv) &
  !$OMP & PRIVATE(ii_ik,ii,ik,ig,igv2,ig2)
  !$OMP DO
  DO ii_ik = 1, 2 * nqbz
     !
     ii = 1 + (ii_ik-1) / nqbz
     ik = 1 +  ii_ik-1  - nqbz*(ii-1)
     !
     DO ig = 1, npw(ik, ii)
        igv2(1:3) = MODULO(igv(1:3,ig,ik,ii), (/dfftt%nr1, dfftt%nr2, dfftt%nr3/))
        ig2 = 1 + igv2(1) + dfftt%nr1 * igv2(2) + dfftt%nr1 * dfftt%nr2 * igv2(3)
        wfc(ig2, 1:nb(ii), 1:npol, ik, ii) = wfc0(ig, 1:nb(ii), 1:npol, ik, ii)
     END DO
     !
  END DO ! ii_ik = 2 * nqbz
  !$OMP END DO
  !$OMP END PARALLEL
  !
  DO ii = 1, 2
     DO ik = 1, nqbz
        DO ipol = 1, npol
           DO ib = 1, nb(ii)
              CALL cfft3d (wfc(1:dfftt%nnr, ib, ipol, ik, ii), &
              & dfftt%nr1, dfftt%nr2, dfftt%nr3, dfftt%nr1, dfftt%nr2, dfftt%nr3, 1, 1)
              !
           END DO
        END DO
     END DO
  END DO
  !
  ! <beta|phi_k>
  !
  IF ( okvan ) THEN
     DO ii = 1, 2
        DO ik = 1, nqbz
           !
           ikv(1) = (ik - 1) / (nq3*nq2)
           ikv(2) = (ik - 1 - ikv(1)*nq2*nq3) / nq3
           ikv(3) =  ik - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
           WHERE(ikv(1:3)*2 + ii - 1 >= (/nq1,nq2,nq3/)) ikv(1:3) = ikv(1:3) - (/nq1,nq2,nq3/)
           !
           xk(1:3) = (REAL(ikv(1:3), dp) + 0.5_dp*REAL(ii-1, dp)) / REAL((/nq1, nq2, nq3/), dp)
           xk(1:3) = MATMUL(bg(1:3, 1:3), xk(1:3))
           !
           DO ig = 1, npw(ik,ii)
              DO igm = 1, ngm
                 IF(ALL(igv(1:3,ig,ik,ii) == mill(1:3,igm))) THEN
                   igk(ig) = igm
                   EXIT
                 END IF
              END DO
           END DO
           CALL init_us_2(npw(ik,ii), igk, xk, vkb)
           DO ipol = 1, npol
              CALL calbec(npw(ik,ii), vkb, wfc0(1:npwx, 1:nb(ii), ipol,ik,ii), &
              &           becwfc(1:nkb, 1:nb(ii), ipol, ik, ii), nb(ii))
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
