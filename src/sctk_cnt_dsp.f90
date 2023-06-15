!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_cnt_dsp
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Compute cnt and dsp
!
SUBROUTINE cnt_and_dsp(n,cnt1,dsp1)
  !
  USE mp_world, ONLY : mpime, nproc
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  INTEGER,INTENT(OUT) :: cnt1, dsp1
  !
  INTEGER :: ii
  INTEGER :: cnt(0:nproc-1), dsp(0:nproc-1)
  !
  cnt(0:nproc-1)        = n / nproc
  cnt(0:mod(n,nproc)-1) = n / nproc + 1
  dsp(0) = 0
  DO ii = 1, nproc - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  END DO
  !
  cnt1 = cnt(mpime)
  dsp1 = dsp(mpime)
  !
END SUBROUTINE cnt_and_dsp
!
! Compute cnt and dsp
!
SUBROUTINE cnt_and_dsp_full(n,cnt,dsp)
  !
  USE mp_world, ONLY : nproc
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  INTEGER,INTENT(OUT) :: cnt(0:nproc - 1), dsp(0:nproc - 1)
  !
  INTEGER :: ii
  !
  cnt(0:nproc-1)        = n / nproc
  cnt(0:mod(n,nproc)-1) = n / nproc + 1
  dsp(0) = 0
  DO ii = 1, nproc - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  END DO
  !
END SUBROUTINE cnt_and_dsp_full
  !
END MODULE sctk_cnt_dsp
