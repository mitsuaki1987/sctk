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

MODULE mp_kel
  !----------------------------------------------------------------------------
  !! Pool groups (processors within a pool of k-points).
  !! Subdivision of image group, used for k-point parallelization.
  !
  USE mp, ONLY : mp_barrier, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  !
  INTEGER :: nproc_band1       = 1
  INTEGER :: nproc_band2  = 1
  INTEGER :: me_band2     = 0
  INTEGER :: me_band1  = 0
  INTEGER :: band1_comm  = 0
  INTEGER :: band2_comm  = 0
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_kel()
    !---------------------------------------------------------------------------
    !! Divide processors (of the "world_comm" group) into "pools".
    !! Requires: \(\text{nproc_band1}\) read from command line,
    !!           \(\text{world_comm}\), typically \(\text{world_comm} =
    !! \text{group}\) of all processors.
    !
    USE mp_world, ONLY : nproc, mpime, world_comm
    IMPLICIT NONE
    !
#if defined (__MPI)
    !
    INTEGER :: nproc_band0, iproc
    !
    nproc_band0 = CEILING(SQRT(DBLE(nproc)))
    !
    DO iproc = nproc_band0, 1, -1
      IF(MOD(nproc, iproc) == 0) THEN
        nproc_band1 = iproc
        nproc_band2 = nproc / iproc
        EXIT
      END IF
    END DO
    !
    ! ... nproc_band1 must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    IF ( nproc_band1 < 1 .OR. nproc_band1 > nproc ) CALL errore( 'mp_start_kel',&
                          'invalid number of pools, out of range', 1 )

    IF ( MOD( nproc, nproc_band1 ) /= 0 ) CALL errore( 'mp_start_kel', &
           'invalid number of pools, nproc /= nproc_band2 * nproc_band1', 1 )
    !
    ! ... me_band1  =  pool index for this processor    ( 0 : nproc_band1 - 1 )
    ! ... me_band2     =  processor index within the pool  ( 0 : nproc_band2 - 1 )
    !
    me_band1 =  mpime / nproc_band2
    me_band2    = MOD( mpime, nproc_band2 )
    !
    CALL mp_barrier( world_comm )
    !
    ! ... the band2_comm communicator is created
    !
    CALL mp_comm_split ( world_comm, me_band1, mpime, band2_comm )
    !
    CALL mp_barrier( world_comm )
    !
    ! ... the band1_comm communicator is created
    !
    CALL mp_comm_split ( world_comm, me_band2, mpime, band1_comm )
    !
#endif
    !
    RETURN
  END SUBROUTINE mp_start_kel
  !
END MODULE mp_kel
