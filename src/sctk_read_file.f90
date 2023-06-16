!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_read_file
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Read elph*.dat
!
SUBROUTINE read_elph()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE modes, ONLY : nmodes
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp, ONLY : mp_bcast, mp_max
  USE disp,  ONLY : nq1, nq2, nq3, x_q, nqs
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE output,        ONLY : fildyn
  USE cell_base, ONLY : bg
  USE io_files, ONLY : prefix, tmp_dir
  !
  USE sctk_val, ONLY : gg0, nqbz, omg0, freq_min, freq_min_ratio
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  INTEGER :: nb0, iq, fi, iqv(3), cnt, dsp, ios
  REAL(dp) :: qvec(3)
  REAL(dp),ALLOCATABLE :: ggc(:,:,:,:,:)
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  CALL start_clock("read_elph")
  !
  ! Read # of k, bands, modes
  !
  IF(ionode) THEN
     fi = find_free_unit()
     OPEN(fi, file = TRIM(tmp_dir) // "_ph0/" // TRIM(prefix) // ".q_1/" &
     &               // TRIM(prefix) // ".elph1", iostat=ios, status="old")
     IF(ios /= 0) CALL errore("read_elph", &
     &                        "Can not open " // TRIM(tmp_dir) // TRIM(prefix) // ".elph1", 1)
     READ(fi,*) nq1, nq2, nq3
     READ(fi,*) elph_nbnd_min, elph_nbnd_max
     READ(fi,*) qvec(1:3)
     READ(fi,*) nmodes
     CLOSE(fi)
  END IF
  !
  CALL mp_bcast(nq1, ionode_id, world_comm)
  CALL mp_bcast(nq2, ionode_id, world_comm)
  CALL mp_bcast(nq3, ionode_id, world_comm)
  CALL mp_bcast(nmodes, ionode_id, world_comm)
  CALL mp_bcast(elph_nbnd_min, ionode_id, world_comm)
  CALL mp_bcast(elph_nbnd_max, ionode_id, world_comm)
  WRITE(stdout,'(7x,"k grid for MEs : ",3(i0,2x))') nq1, nq2, nq3
  WRITE(stdout,'(7x,"Number of modes : ",i0)') nmodes
  !
  ! Search irreducible k points
  !
  fildyn       = 'matdyn'
  CALL q_points()
  nqbz = nq1 * nq2 * nq3
  !
  CALL cnt_and_dsp(nqs,cnt,dsp)
  !
  ALLOCATE(gg0(nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,nqbz,cnt), &
  &        omg0(nmodes,cnt))
  !
  ! Read omega & gg
  !
  ALLOCATE(ggc(2,nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,nqbz))
  !
  DO iq = 1, cnt
     !
     fi = find_free_unit()
     OPEN(fi, file = TRIM(tmp_dir) // "_ph0/" // TRIM(prefix) // ".q_" // TRIM(int_to_char(dsp+iq)) // "/" &
     &               // TRIM(prefix) // ".elph" // TRIM(int_to_char(dsp+iq)), &
     &    iostat=ios, status="old") !, form = 'unformatted')
     IF(ios /= 0) CALL errore("read_elph", &
     &  "Can not open " // TRIM(tmp_dir) // TRIM(prefix) // ".elph"//TRIM(int_to_char(dsp+iq)), dsp+iq)
     !
     READ(fi,*) iqv(1:3)
     IF(.NOT. ALL(iqv(1:3) == (/nq1, nq2, nq3/))) &
     & CALL errore("read_elph", "k grid. iq = ", dsp + iq)
     !
     READ(fi,*) iqv(1:2)
     IF(iqv(1) /= elph_nbnd_min) &
     & CALL errore("read_elph", "First band. iq = ", dsp + iq)
     !
     IF(iqv(2) /= elph_nbnd_max) &
     & CALL errore("read_elph", "Last band. iq = ", dsp + iq)
     !
     READ(fi,*) qvec(1:3)
     qvec(1:3) = MATMUL(bg(1:3,1:3), qvec(1:3))
     !
     IF(ANY(ABS(qvec(1:3) - x_q(1:3,iq + dsp)) > 1e-5_dp)) &
     & CALL errore("read_elph", "qvec. iq = ", dsp + iq)
     !
     READ(fi,*) nb0
     IF(nb0 /= nmodes) &
     & CALL errore("read_elph", "# of modes. iq = ", dsp + iq)
     !
     READ(fi,*) omg0(1:nmodes,iq)
     READ(fi,*) ggc(1:2,1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz)
     gg0(               1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz,iq) &
     & = ggc(         1,1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz)**2 &
     & + ggc(         2,1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz)**2
     !
     CLOSE(fi)
     !
  END DO
  !
  DEALLOCATE(ggc)
  !
  IF(freq_min_ratio > 0.0) THEN
     freq_min = freq_min_ratio * MAXVAL(omg0(1:nmodes,1:cnt))
     CALL mp_max(freq_min, world_comm)
  END IF
  !
  CALL stop_clock("read_elph")
  !
END SUBROUTINE read_elph
!
! Read Screened Coulomb matrix elements
!
SUBROUTINE read_Coulomb()
  !
  USE wvfct, ONLY : nbnd
  USE kinds, ONLY : DP
  USE mp, ONLY : mp_bcast
  USE disp,  ONLY : nq1, nq2, nq3, x_q, nqs
  USE cell_base, ONLY : bg
  USE mp_world, ONLY : world_comm, mpime
  USE io_files, ONLY : prefix, tmp_dir
  !
  USE sctk_val, ONLY : nqbz, Vc0, lsf, nci, lz_coulomb
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  INTEGER :: fi, iq, iqv(3), nb0, cnt, dsp, ios
  REAL(dp) :: qvec(3)
  !
  CHARACTER(LEN=6) :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  CALL start_clock("read_Coulomb")
  !
  IF(mpime==0) THEN
     fi = find_free_unit()
     OPEN(fi, file = TRIM(tmp_dir) // TRIM(prefix) // ".vel1", &
     &    form = 'unformatted', iostat=ios, status="old")
     IF(ios /= 0) CALL errore("read_Coulomb", &
     &                        "Can not open " // TRIM(tmp_dir) // TRIM(prefix) // ".vel1", 1)
     READ(fi) iqv(1:3)
     READ(fi) nb0
     READ(fi) qvec(1:3)
     READ(fi) nci
     WRITE(*,'(7x,"Dimention of Chebyshev interpolation : ",i0)') nci
     CLOSE(fi)
  END IF
  !
  CALL mp_bcast(nci, 0, world_comm)
  !
  CALL cnt_and_dsp(nqs,cnt,dsp)
  ALLOCATE(vc0(nci*lsf,nbnd,nbnd,nqbz,cnt))
  !
  DO iq = 1, cnt
     !
     fi = find_free_unit()
     OPEN(fi, file = TRIM(tmp_dir) // TRIM(prefix) // ".vel"//TRIM(int_to_char(iq+dsp)), &
     &    form = 'unformatted', iostat=ios, status="old")
     IF(ios /= 0) CALL errore("read_Coulomb", &
     &  "Can not open " // TRIM(tmp_dir) // TRIM(prefix) // ".vel"//TRIM(int_to_char(iq+dsp)), iq+dsp)
     !
     READ(fi) iqv(1:3)
     IF(.NOT. ALL(iqv(1:3) == (/nq1,nq2,nq3/))) &
     &  CALL errore("read_coulomb", "kgrid. iq = ", dsp + iq)
     !
     READ(fi) nb0
     IF(nbnd /= nb0) &
     &  CALL errore("read_coulomb", "# of bands. iq = ", dsp + iq)
     !
     READ(fi) qvec(1:3)
     qvec(1:3) = MATMUL(bg(1:3,1:3), qvec(1:3))
     IF(ANY(ABS(qvec(1:3) - x_q(1:3,iq + dsp)) > 1e-5_dp)) &
     &  CALL errore("read_coulomb", "q point. iq = ", dsp + iq)
     !
     READ(fi) nb0
     IF(nci /= nb0) &
     &  CALL errore("read_coulomb", "Chebyshev interpolation. iq = ", dsp + iq)
     !
     READ(fi) vc0(1:nci,1:nbnd,1:nbnd,1:nqbz,iq)
     !
     IF(lsf==2) THEN
        !
        READ(fi) vc0(nci+1:nci*2,1:nbnd,1:nbnd,1:nqbz,iq)
        !
        ! V1 = V_C + V_S
        !
        vc0(        1:nci,  1:nbnd,1:nbnd,1:nqbz,iq) &
        & = vc0(    1:nci,  1:nbnd,1:nbnd,1:nqbz,iq) &
        & + vc0(nci+1:nci*2,1:nbnd,1:nbnd,1:nqbz,iq)
        !
        IF(lz_coulomb) THEN
           !
           ! V2 = V_S - V_C = 2 V_S - (V_C + V_S))
           !
           vc0(    nci+1:nci*2,1:nbnd,1:nbnd,1:nqbz,iq) &
           & = vc0(nci+1:nci*2,1:nbnd,1:nbnd,1:nqbz,iq) * 2.0_dp &
           & - vc0(    1:nci  ,1:nbnd,1:nbnd,1:nqbz,iq)
           !
        END IF
        !
     END IF
     !
     CLOSE(fi)
     !
  END DO
  !
  CALL stop_clock("read_Coulomb")
  !
END SUBROUTINE read_Coulomb
!
! Read from data-file.xml
!
SUBROUTINE read_a2fsave()
  !
  USE kinds, ONLY : DP
  USE parameters, ONLY : npk
  USE mp_world, ONLY : world_comm
  USE io_files, ONLY : prefix, tmp_dir
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp, ONLY : mp_bcast
  USE start_k, ONLY : nk1, nk2, nk3
  USE fermisurfer_common, ONLY : rotate_k_fs
  USE wvfct, ONLY : nbnd, wg, et
  USE klist, ONLY : nks, xk, wk, nelec, nkstot
  USE ener, ONLY : ef
  USE cell_base, ONLY : at, bg
  USE lsda_mod, ONLY : nspin, isk
  USE ktetra, ONLY : opt_tetra_init, opt_tetra_weights, tetra_type
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev
  !
  INTEGER :: fi, ik, i1, i2, i3, s_dummy(3,3,48), t_rev_dummy(48)
  INTEGER,ALLOCATABLE :: equiv(:,:,:)
  REAL(dp),ALLOCATABLE :: et0(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  !
  IF ( ionode ) THEN
     fi = find_free_unit()
     OPEN(fi, file = TRIM(tmp_dir) // TRIM(prefix) // ".a2Fsave")
     READ(fi,*) nbnd, nks
     ALLOCATE(et0(nbnd,nks))
     READ(fi,*) et0(1:nbnd,1:nks)
     READ(fi,*) xk(   1:3,1:nks)
     READ(fi,*) wk(       1:nks)
     READ(fi,*) nk1, nk2, nk3
     CLOSE(fi)
  END IF
  !
  CALL mp_bcast (nbnd, ionode_id, world_comm)
  CALL mp_bcast (nk1, ionode_id, world_comm)
  CALL mp_bcast (nk2, ionode_id, world_comm)
  CALL mp_bcast (nk3, ionode_id, world_comm)
  CALL mp_bcast (nks, ionode_id, world_comm)
  IF(.NOT. ionode) ALLOCATE(et0(nbnd,nks))
  CALL mp_bcast (et0, ionode_id, world_comm)
  CALL mp_bcast (xk, ionode_id, world_comm)
  CALL mp_bcast (wk, ionode_id, world_comm)
  !
  ! ... Find equivalent k point in irr-BZ for whole BZ
  !
  WRITE(stdout,'(5x,"Dense k-grid : ",3(i0,2x))') nk1, nk2, nk3
  WRITE(stdout,'(5x,"Number of k(dense) : ",i0)') nks
  ALLOCATE(equiv(nk1, nk2, nk3))
  !
  IF(ALLOCATED(et)) DEALLOCATE(et)
  IF(ALLOCATED(wg)) DEALLOCATE(wg)
  !
  ALLOCATE(et(nbnd,nks), wg(nbnd,nks))
  et(1:nbnd, 1:nks) = et0(1:nbnd, 1:nks)
  tetra_type = 2
  CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, 0, 0, 0, &
  &                   nk1, nk2, nk3, nks, xk, 1)
  !
  CALL opt_tetra_weights (nks, nspin, nbnd, nelec, et, ef, wg, 0, isk)
  CALL rotate_k_fs(equiv)
  !
  ! ... Map e_k into whole BZ (Measured from E_F)
  !
  DEALLOCATE(et, wg)
  nks = nk1 * nk2 * nk3
  nkstot = nks * 2
  ALLOCATE(et(nbnd,nks*2), wg(nbnd,nks))
  ik = 0
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        DO i3 = 1, nk3
           ik = ik + 1
           et(1:nbnd,ik) = et0(1:nbnd, equiv(i1,i2,i3))
           xk(1:3,ik) = REAL((/i1, i2, i3/), DP) / REAL((/nk1, nk2, nk3/), DP)
           WHERE((/i1, i2, i3/)*2 >= (/nk1, nk2, nk3/)) xk(1:3,ik) = xk(1:3,ik) - 1.0_dp
           xk(1:3,ik) = MATMUL(bg(1:3,1:3), xk(1:3,ik))
           wk(ik) = 1.0_dp / REAL(nks, DP)
        END DO
     END DO
  END DO
  !
  ! ... Find Fermi energy
  !
  s_dummy(1:3,1:3,1:48) = 0
  t_rev_dummy(1:48) = 0
  DO i1 = 1, 3
     s_dummy(i1,i1,1:48) = 1
  END DO
  CALL opt_tetra_init(1, s_dummy, .False., t_rev_dummy, at, bg, npk, &
  &                   0, 0, 0, nk1, nk2, nk3, nks, xk, 1)
  !
  DEALLOCATE(et0,equiv)
  !
END SUBROUTINE read_a2fsave
!
SUBROUTINE degenerated_band()
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : nbnd, et
  USE sctk_val, ONLY : degen, ndegen
  USE klist, ONLY : nks
  !
  INTEGER :: ibnd, iks, ii, nksb2
  REAL(DP) :: et0, et1
  !
  nksb2 = nks / 2
  !
  ALLOCATE(ndegen(nksb2,2), degen(2,nbnd,nksb2,2))
  !
  DO ii = 1, 2
     DO iks = 1, nksb2
        !
        ndegen(iks,ii) = 1
        degen(1,ndegen(iks,ii),iks,ii) = 1
        et0 = et(1,iks + nksb2*(ii-1))
        !
        DO ibnd = 2, nbnd
           !
           et1 = et(ibnd,iks + nksb2*(ii-1))
           !
           IF(ABS(et0 - et1) > 1.0e-7_dp) THEN
              degen(2,ndegen(iks,ii),iks,ii) = ibnd - 1
              et0 = et1
              ndegen(iks,ii) = ndegen(iks,ii) + 1
              degen(1,ndegen(iks,ii),iks,ii) = ibnd
           END IF
           !
        END DO
        !
        degen(2,ndegen(iks,ii),iks,ii) = nbnd
        !
     END DO
  END DO
  !
END SUBROUTINE degenerated_band
!
END MODULE sctk_read_file
