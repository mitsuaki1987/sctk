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
  USE mp_pools, ONLY : intra_pool_comm
  USE mp, ONLY : mp_max, mp_sum
  USE disp, ONLY : nq1, nq2, nq3
  USE io_global, ONLY : stdout
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : nqbz, wfc, igv, npw, nkpe
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  INTEGER :: fi, ik, ii, ik_g, ngw0
  INTEGER :: cnt, dsp
  CHARACTER(LEN=320) :: filename
  !
  INTEGER :: ik0, ispin, npol0, ierr, ikv(3), ibnd, npwxtot
  REAL(dp) :: scalef, xk0(3), bvec0(3,3), xk1(3)
  LOGICAL :: gamma_only
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  fi = find_free_unit()
  !
  CALL cnt_and_dsp(nqbz,cnt,dsp)
  nkpe = cnt
  CALL mp_max(nkpe, world_comm)
  WRITE(stdout,'(7x,"k-point per PE : ",i0)') nkpe
  !
  npwxtot = npwx
  CALL mp_sum(npwxtot, intra_pool_comm)
  !
  WRITE(stdout,'(7x,"Total RAM for input WFC : ",e10.2," GB")') &
  &     REAL(npwxtot,dp)*REAL(npol*2,dp)*REAL(nbnd,dp)*REAL(nqbz,dp)*16.0e-9_dp
  ALLOCATE(wfc(npwxtot,npol,nbnd,nkpe,2), igv(3,npwxtot,nkpe,2), npw(nkpe,2))
  igv(1:3,1:npwxtot,            1:nkpe,1:2) = 0
  wfc(    1:npwxtot,1:npol,1:nbnd,1:nkpe,1:2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Read Wfc(k)
  !
  DO ii = 1, 2
     !
     DO ik = 1, cnt
        !
        ik_g = nqbz * (ii - 1) + dsp + ik
        filename = TRIM(tmp_dir) // TRIM(prefix) // '.save/' // 'wfc' // TRIM(int_to_char(ik_g))
        !
        OPEN(unit=fi, file=trim(filename)//'.dat', form='unformatted', status='old', iostat=ierr)
        READ(fi) ik0, xk0(1:3), ispin, gamma_only, scalef
        READ(fi) ngw0, npw(ik,ii), npol0, ibnd
        READ(fi) bvec0(1:3,1:3)
        READ(fi) igv(1:3,1:npw(ik,ii),ik,ii)
        DO ibnd = 1, nbnd
           READ(fi) wfc(1:npw(ik,ii),1:npol,ibnd,ik,ii)
        END DO
        CLOSE(unit=fi, status='keep')
        !
        ! Verify k-point
        !
        ik_g = dsp + ik
        ikv(1) = (ik_g - 1) / (nq3*nq2)
        ikv(2) = (ik_g - 1 - ikv(1)*nq2*nq3) / nq3
        ikv(3) =  ik_g - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
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
     END DO ! ik
     !
  END DO ! ii
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
  USE fft_scalar, ONLY : cfft3d
  USE fft_support, only: good_fft_order
  USE mp_world, ONLY : world_comm
  USE mp_pools, ONLY : intra_pool_comm
  USE io_global, ONLY : stdout
  USE noncollin_module, ONLY : npol
  !
  USE sctk_val, ONLY : igv, nf, nftot, npw, wfc, &
  &                    wfc1, wfc2, wfc1q, wfc2q, nqbz, nkpe
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ik, ib, igv2(3), ig2, cnt, dsp, ii, igmin, igmax, npwxtot, ipol
  !
  npwxtot = npwx
  CALL mp_sum(npwxtot, intra_pool_comm)
  CALL cnt_and_dsp(nqbz,cnt,dsp)
  !
  ! FFT grid
  !
  DO ii = 1, 3
     igmin = MINVAL(igv(1:3, 1:npwxtot, 1:cnt, 1:2))
     igmax = MAXVAL(igv(1:3, 1:npwxtot, 1:cnt, 1:2))
     CALL mp_max(igmin, world_comm)
     CALL mp_max(igmax, world_comm)     
     nf(ii) = good_fft_order(igmax - igmin + 1)
  END DO
  nftot = PRODUCT(nf(1:3))
  WRITE(stdout,'(7x,"Grid :",3(2x,i0))') nf(1:3)
  !
  WRITE(stdout,'(7x,"Total RAM for WFC(r) : ",e10.2," GB")') &
  &     REAL(nftot,dp)*REAL(npol*4,dp)*REAL(nbnd,dp)*REAL(nqbz,dp)*16.0e-9_dp
  !
  ALLOCATE(wfc1(nftot, npol, nbnd, nkpe), wfc2(nftot, npol, nbnd, nkpe))
  ALLOCATE(wfc1q(nftot, npol, nbnd, nkpe), wfc2q(nftot, npol, nbnd, nkpe))
  wfc1(1:nftot,1:npol,1:nbnd,1:nkpe) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  wfc2(1:nftot,1:npol,1:nbnd,1:nkpe) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
  !
  ! Fourier trans. wfc(G) -> wfc(r)
  ! 
  DO ik = 1, cnt
     !
     DO ig = 1, npw(ik,1)
        igv2(1:3) = MODULO(igv(1:3,ig,ik,1), nf(1:3))
        ig2 = 1 + igv2(1) + nf(1) * igv2(2) + nf(1) * nf(2) * igv2(3)
        wfc1(ig2,1:npol, 1:nbnd, ik) = wfc(ig, 1:npol, 1:nbnd, ik, 1)
     END DO
     !
     DO ig = 1, npw(ik,2)
        igv2(1:3) = MODULO(igv(1:3,ig,ik,2), nf(1:3))
        ig2 = 1 + igv2(1) + nf(1) * igv2(2) + nf(1) * nf(2) * igv2(3)
        wfc2(ig2, 1:npol, 1:nbnd, ik) = wfc(ig, 1:npol, 1:nbnd, ik, 2)
     END DO
     !
     DO ib = 1, nbnd
        DO ipol = 1, npol
           !
           CALL cfft3d (wfc1(1:nftot, ipol, ib, ik), &
           & nf(1), nf(2), nf(3), nf(1), nf(2), nf(3), 1, 1)
           !
           CALL cfft3d(wfc2(1:nftot, ipol, ib, ik), &
           & nf(1), nf(2), nf(3), nf(1), nf(2), nf(3), 1, 1)
           !
        END DO
     END DO
     !
  END DO
  !
  DEALLOCATE(wfc, igv, npw)
  !
END SUBROUTINE fft_wfc
!
END MODULE sctk_wfc
