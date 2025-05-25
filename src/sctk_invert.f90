!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_invert
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Inversion of matrices
!
SUBROUTINE invert(matrix)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : nproc, mpime, world_comm
  USE mp, ONLY : mp_sum
  USE sctk_val, ONLY : ngv0, ngv1, ngv, nmf
  !
  IMPLICIT NONE
  !
  COMPLEX(dp),INTENT(INOUT) :: matrix(ngv,ngv0:ngv1,0:nmf)
  !
  INTEGER :: imf, ig, jg, imaxpiv(nproc)
  REAL(dp) :: maxpiv(nproc)
  COMPLEX(dp) :: key1(ngv), key2(ngv), piv
  COMPLEX(dp),ALLOCATABLE :: matrix1(:,:), matrix2(:,:)
  !
  CALL start_clock("invert")
  !
  ALLOCATE(matrix1(ngv,ngv0:ngv1), matrix2(ngv,ngv0:ngv1))
  !
  DO imf = 0, nmf
     !
     matrix1(1:ngv, ngv0:ngv1) = matrix(1:ngv, ngv0:ngv1,imf)
     !
     matrix2(1:ngv,ngv0:ngv1) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
     DO ig = ngv0, ngv1
        matrix2(ig,ig) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
     END DO
     !
     DO ig = 1, ngv
        !
        ! Partial pivoting
        !
        jg = MAX(ig, ngv0)
        maxpiv( 1:nproc) = 0.0_dp
        imaxpiv(1:nproc) = 0
        IF(jg > ngv1) THEN
           maxpiv( mpime+1) = - 1e10_dp
           imaxpiv(mpime+1) =   1
        ELSE
           maxpiv( mpime+1) = MAXVAL(REAL(CONJG(matrix1(ig, jg:ngv1)) &
           &                                  * matrix1(ig, jg:ngv1), dp ))
           imaxpiv(mpime+1) = MAXLOC(REAL(CONJG(matrix1(ig, jg:ngv1)) &
           &                                  * matrix1(ig, jg:ngv1), dp ), 1) &
           &                + jg - 1
        END IF
        !
        CALL mp_sum(maxpiv, world_comm)
        CALL mp_sum(imaxpiv, world_comm)
        jg = imaxpiv(MAXLOC(maxpiv, 1))
        !
        ! Exchange matrix1
        !
        key1(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        key2(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        !
        IF(ngv0 <= ig  .AND. ig  <= ngv1) key1(1:ngv) = matrix1(1:ngv,ig)
        IF(ngv0 <= jg  .AND. jg  <= ngv1) key2(1:ngv) = matrix1(1:ngv,jg)
        !
        CALL mp_sum( key1, world_comm )
        !
        CALL mp_sum( key2, world_comm )
        !
        IF(ngv0 <= ig  .AND. ig  <= ngv1) matrix1(1:ngv,ig) = key2(1:ngv)
        IF(ngv0 <= jg  .AND. jg  <= ngv1) matrix1(1:ngv,jg) = key1(1:ngv)
        !
        ! Exchange matrix2
        !
        key1(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        key2(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        !
        IF(ngv0 <= ig  .AND. ig  <= ngv1) key1(1:ngv) = matrix2(1:ngv,ig)
        IF(ngv0 <= jg  .AND. jg  <= ngv1) key2(1:ngv) = matrix2(1:ngv,jg)
        !
        CALL mp_sum( key1, world_comm )
        !
        CALL mp_sum( key2, world_comm )
        !
        IF(ngv0 <= ig  .AND. ig  <= ngv1) matrix2(1:ngv,ig) = key2(1:ngv)
        IF(ngv0 <= jg  .AND. jg  <= ngv1) matrix2(1:ngv,jg) = key1(1:ngv)
        !
        ! Ordinary Gauss-Jordan
        !
        IF(ngv0 <= ig .AND. ig <= ngv1) THEN
           !
           piv = CMPLX(1.0_dp, 0.0_dp, KIND=dp) / matrix1(ig, ig)
           !
           IF(ABS(piv) < 1e-12_dp) THEN
              CALL errore ('invert', 'Singular imf = ', imf)
           END IF
           !
           CALL zscal(ngv, piv, matrix1(1:ngv, ig), 1)
           CALL zscal(ngv, piv, matrix2(1:ngv, ig), 1)
           CALL zcopy(ngv, matrix1(1:ngv, ig), 1, key1, 1)
           CALL zcopy(ngv, matrix2(1:ngv, ig), 1, key2, 1)
           !  
        ELSE
           !
           key1(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
           key2(1:ngv) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
           !
        END IF ! IF(ngv0 <= ig .AND. ig <= ngv1) 
        !
        CALL mp_sum( key1, world_comm )
        !
        CALL mp_sum( key2, world_comm )
        !
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP & SHARED(ngv0,ngv1,ig,matrix1,matrix2,ngv,key1,key2) &
        !$OMP & PRIVATE(jg,piv)
        !
        !$OMP DO
        DO jg = ngv0, ngv1
           !
           IF(jg == ig) CYCLE
           !
           piv = - matrix1(ig, jg)
           CALL zaxpy(ngv, piv, key1, 1, matrix1(1:ngv, jg), 1)
           CALL zaxpy(ngv, piv, key2, 1, matrix2(1:ngv, jg), 1)
           !
        END DO ! jg
        !$OMP END DO
        !$OMP END PARALLEL
        !
     END DO ! ig
     !
     matrix(1:ngv,ngv0:ngv1, imf) = matrix2(1:ngv,ngv0:ngv1)
     !
  END DO ! imf
  !
  DEALLOCATE(matrix1,matrix2)
  !
  CALL stop_clock("invert")
  !
END SUBROUTINE invert
!
! Invertion of matrices
!
SUBROUTINE hermite()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm, nproc
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp_full
  USE noncollin_module, ONLY : npol
  USE parallel_include
  !
  USE sctk_val, ONLY : nmf, wscr, ngv, ngv0, ngv1
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, jg, imf, ierr, ip, ipol, &
  &          cnt(0:nproc-1), dsp(0:nproc-1), cnt0, dsp0
  COMPLEX(DP),ALLOCATABLE :: wsnd(:,:), wrcv(:)
  !
  ! Hermite conjugate
  !
  ALLOCATE(wsnd(ngv0:ngv1,ngv), wrcv((ngv1-ngv0+1)*ngv))
  CALL cnt_and_dsp_full(ngv,cnt,dsp)
  cnt(0:nproc-1) = cnt(0:nproc-1) * (ngv1 - ngv0 + 1)
  dsp(0:nproc-1) = dsp(0:nproc-1) * (ngv1 - ngv0 + 1)
  !
  DO ipol = 2, 2*npol
     DO imf = 0, nmf
        !
        DO ig = ngv0, ngv1
           wsnd(ig,1:ngv) = CONJG(wscr(1:ngv,ig,imf,ipol))
        END DO
        !
#if defined (__MPI)
        CALL MPI_ALLTOALLV(wsnd, cnt, dsp, MPI_DOUBLE_COMPLEX, &
        &                  wrcv, cnt, dsp, MPI_DOUBLE_COMPLEX, world_comm, ierr)
        IF(ierr /= 0) CALL errore('lambda_sf', 'mpi_alltoallv', ierr)
        !
        DO ip = 0, nproc - 1
           !
           cnt0 = cnt(ip) / (ngv1 - ngv0 + 1)
           dsp0 = dsp(ip) / (ngv1 - ngv0 + 1)
           !
           !$OMP PARALLEL DEFAULT(NONE) &
           !$OMP & SHARED(ngv0,ngv1,dsp0,cnt0,wscr,wrcv,imf,ipol,dsp,ip) &
           !$OMP & PRIVATE(ig,jg)
           !$OMP DO
           DO ig = ngv0, ngv1
              jg = dsp(ip) + (ig - ngv0) * cnt0
              wscr(   dsp0+1:dsp0+cnt0, ig, imf, ipol) &
              & = wrcv(jg+1:  jg+cnt0)
           END DO
           !$OMP END DO
           !$OMP END PARALLEL
           !
        END DO
#else
        wscr(1:ngv,1:ngv,imf,ipol) = wsnd(1:ngv,1:ngv)
#endif
        !
     END DO ! imf = 0, nmf
  END DO !ipol = 1, 2*npol - 1
  !
  DEALLOCATE(wsnd, wrcv)
  !
END SUBROUTINE hermite
!
! Compute maximum-absolute eigenvalue with power method
!
SUBROUTINE eigmax_power(matrix, eigval, iter)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime, world_comm
  USE mp, ONLY : mp_bcast, mp_sum
  USE sctk_val, ONLY : ngv0, ngv1, ngv
  !
  COMPLEX(dp),INTENT(IN) :: matrix(ngv,ngv0:ngv1)
  REAL(dp),INTENT(OUT) :: eigval
  !
  INTEGER :: iter
  REAL(dp) :: norm, res
  !
  REAL(dp),ALLOCATABLE :: rvec(:)
  COMPLEX(dp),ALLOCATABLE :: vec_i(:), vec_o(:)
  !
  ALLOCATE(vec_i(ngv), vec_o(ngv), rvec(ngv))
  !
  ! Initial guess is random vector
  !
  CALL random_seed()
  CALL random_number(rvec)
  CALL mp_bcast(rvec, 0, world_comm)
  norm  = REAL(DOT_PRODUCT(rvec(1:ngv), rvec(1:ngv)), dp)
  vec_i(1:ngv) = rvec(1:ngv) / SQRT(norm)
  !
  DO iter = 1, 500
    !
    call ZGEMV("N", ngv, ngv1-ngv0+1, 1.0_dp, matrix, ngv, vec_i(ngv0:ngv1), 1, 0.0_dp, vec_o(   1:ngv ), 1)
    !
    CALL mp_sum(vec_o, world_comm)
    !
    res   = REAL(DOT_PRODUCT(vec_i(1:ngv), vec_o(1:ngv)), dp) - eigval
    norm  = REAL(DOT_PRODUCT(vec_o(1:ngv), vec_o(1:ngv)), dp)
    eigval = eigval + res
    !
    vec_i(1:ngv) = vec_o(1:ngv) / SQRT(norm)
    !
    IF(ABS(res) < 1.0e-5_dp .AND. iter > 1) EXIT
    !
  END DO
  !
  DEALLOCATE(vec_i, vec_o, rvec)
  !
END SUBROUTINE
!
END MODULE sctk_invert
