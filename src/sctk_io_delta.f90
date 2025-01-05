!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_io_delta
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Define initial delta
!
SUBROUTINE read_delta()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime, world_comm
  USE mp, ONLY : mp_bcast
  USE io_global, ONLY : ionode_id
  USE io_files, ONLY : prefix, tmp_dir
  USE sctk_val, ONLY : bindx, delta, dk, kindx, ngap, ngap1, ngap2, xi
  !
  IMPLICIT NONE
  !
  INTEGER :: it, fi = 10, is
  REAL(dp) :: Z0
  CHARACTER(1) :: tmp
  !
  IF(mpime == 0) THEN
     !
     OPEN(fi, file = TRIM(tmp_dir) // TRIM(prefix) // '.scgap',status="old", action = 'read',iostat = is)
     IF(is /= 0) CALL errore("read_delta", &
     &                       "Can not open" // TRIM(tmp_dir) // TRIM(prefix) // ".scgap", 1)
     !
     READ(fi,*) tmp, ngap1, ngap2
     !
     ngap = MAX(ngap1, ngap2)
     !
     WRITE(*,'(7x,"Number of total points for gap equation : ",2(i0,2x))') ngap1, ngap2
     ALLOCATE(xi(ngap,2), delta(ngap,2), dk(ngap,2), kindx(ngap,2), bindx(ngap,2))
     !
     DO it = 1, ngap1
        !
        READ(fi,*) xi(it,1), delta(it,1), Z0, dk(it,1), kindx(it,1), bindx(it,1)
        !
     END DO
     !
     DO it = 1, ngap2
        !
        READ(fi,*) xi(it,2), delta(it,2), Z0, dk(it,2), kindx(it,2), bindx(it,2)
        !
     END DO
     !
     CLOSE(fi)
     !
  END IF
  !
  CALL mp_bcast(ngap,  ionode_id, world_comm )
  CALL mp_bcast(ngap1, ionode_id, world_comm )
  CALL mp_bcast(ngap2, ionode_id, world_comm )
  IF(mpime /= 0) ALLOCATE(xi(ngap,2), delta(ngap,2), dk(ngap,2), kindx(ngap,2), bindx(ngap,2))
  CALL mp_bcast(xi,    ionode_id, world_comm )
  CALL mp_bcast(delta, ionode_id, world_comm )
  CALL mp_bcast(dk,    ionode_id, world_comm )
  CALL mp_bcast(kindx, ionode_id, world_comm )
  CALL mp_bcast(bindx, ionode_id, world_comm )
  !
END SUBROUTINE read_delta
!
! Output to file (delta)
!
SUBROUTINE out_delta(fname)
  !
  USE mp_world, ONLY : mpime
  USE sctk_val, ONLY : bindx, delta, dk, kindx, ngap1, ngap2, xi, Z
  IMPLICIT NONE
  !
  CHARACTER(*),INTENT(IN) :: fname
  !
  INTEGER :: it, fo = 20
  !
  IF(mpime == 0) THEN
    !
    OPEN(fo, file = fname)
    !
    WRITE(fo,*) "#", ngap1, ngap2
    !
    WRITE(fo,*) ""
    !
    DO it = 1, ngap1
      WRITE(fo,'(4e25.15,2i8)') xi(it,1), delta(it,1), Z(it,1), dk(it,1), &
      &                       kindx(it,1), bindx(it,1)
    END DO
    !
    WRITE(fo,*) ""
    !
    DO it = 1, ngap2
      WRITE(fo,'(4e25.15,2i8)') xi(it,2), delta(it,2), Z(it,2), dk(it,2), &
      &                         kindx(it,2), bindx(it,2)
    END DO
    !
    CLOSE(fo)
    !
  END IF
  !
END SUBROUTINE out_delta
!
! write QPDOS
!
SUBROUTINE write_dos()
  !
  USE kinds, ONLY : DP
  USE sctk_val, ONLY : e0, ne, sdos
  USE constants, ONLY : RYTOEV
  USE io_files, ONLY : prefix, tmp_dir
  !
  INTEGER :: ie, fo = 20
  !
  OPEN(fo, file = TRIM(tmp_dir) // TRIM(prefix) // '.qpdos')
  !
  DO ie = 1, ne
     WRITE(fo,*) e0(ie) * RYTOEV * 1.0e3_dp, sdos(ie)
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE write_dos
!
! Output .frmsf file
!
SUBROUTINE output_frmsf(eig,mat,filename)
  !
  USE kinds,     ONLY : DP
  USE start_k,   ONLY : nk1, nk2, nk3
  USE wvfct, ONLY : nbnd
  USE ener, ONLY : ef
  USE fermisurfer_common,   ONLY : b_low, b_high, write_fermisurfer
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: eig(nbnd,        nk3,nk2,nk1), &
  &                      mat(b_low:b_high,nk3,nk2,nk1)
  CHARACTER(*),INTENT(IN) :: filename
  !
  INTEGER :: i1, i2
  REAL(dp) :: eig1(b_low:b_high,nk1,nk2,nk3), mat1(b_low:b_high,nk1,nk2,nk3)
  !
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        eig1(b_low:b_high,i1,i2,1:nk3) = eig(b_low:b_high,1:nk3,i2,i1) - ef
        mat1(b_low:b_high,i1,i2,1:nk3) = mat(b_low:b_high,1:nk3,i2,i1)
     END DO
  END DO
  !
  CALL write_fermisurfer(eig1, mat1, filename)  
  !
END SUBROUTINE output_frmsf
!
END MODULE sctk_io_delta
