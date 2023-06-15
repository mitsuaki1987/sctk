!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_stdin
  !
  IMPLICIT NONE
  !
CONTAINS
!>
!> Read "control" namelist
!>
SUBROUTINE stdin_control()
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_bcast
  USE io_files, ONLY : prefix, tmp_dir
  USE input_parameters, ONLY : calculation
  !
  IMPLICIT NONE
  !
  CHARACTER(256) :: outdir
  CHARACTER(256), EXTERNAL :: trimcheck
  NAMELIST /control/ prefix, outdir, calculation
  !
  prefix = 'pwscf'
  CALL get_environment_variable('ESPRESSO_TMPDIR', outdir)
  IF(TRIM(outdir) == ' ') outdir = './'
  calculation = 'kel'
  !
  IF(ionode) THEN
     !
     CALL input_from_file ( )
     READ(5,control,err=100)
     tmp_dir = trimcheck(outdir)
     !
  END IF
  !
  CALL mp_bcast(prefix,  ionode_id, world_comm)
  CALL mp_bcast(tmp_dir, ionode_id, world_comm)
  CALL mp_bcast(calculation, ionode_id, world_comm)
  !
  IF(TRIM(calculation) /= "kel" .AND. &
  &  TRIM(calculation) /= "scdft" .AND. &
  &  TRIM(calculation) /= "scdft_tc" .AND. &
  &  TRIM(calculation) /= "qpdos" .AND. &
  &  TRIM(calculation) /= "deltaf" .AND. &
  &  TRIM(calculation) /= "ultrasonic" .AND. &
  &  TRIM(calculation) /= "lambda_mu_k") THEN
     WRITE(*,*) "calculation = ", TRIM(calculation)
     CALL errore ('stdin_control', 'Invalid input for keyword calculation.', 1)
  END IF
  !
  RETURN
  !
100 CALL errore ('stdin_control', 'reading namelist CONTROL', 1)
  !
END SUBROUTINE stdin_control
!>
!> Read from STDIN for K_el
!>
SUBROUTINE stdin_Kel()
  !
  USE io_global, ONLY : ionode, ionode_id
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_bcast
  USE gvecw, ONLY : ecutwfc
  USE control_ph, ONLY : start_q, last_q
  USE start_k, ONLY : nk1, nk2, nk3
  USE disp,  ONLY : nq1, nq2, nq3, nqs
  USE output,        ONLY : fildyn
  !
  USE sctk_val, ONLY : laddxc, nqbz, lsf, nci
  !
  IMPLICIT NONE
  !
  NAMELIST /kel/ start_q, last_q, nci, laddxc, ecutwfc, nq1, nq2, nq3, lsf
  !
  IF(ionode) THEN
     !
     nci = 5
     start_q = 1
     last_q = 0
     laddxc = 0
     nq1 = nk1
     nq2 = nk2
     nq3 = nk3
     lsf = 0
     !
     READ(5,kel,err=100)
     !
     WRITE(*,'(7x,"               q grid : ",3(i0,2x))') nq1, nq2, nq3
     WRITE(*,'(7x,"  # of Chebyshev int. : ",i0)') nci
     WRITE(*,'(7x,"          The first q : ",i0)') start_q
     WRITE(*,'(7x,"           The last q : ",i0)') last_q
     WRITE(*,'(7x,"               laddxc : ",i0)') laddxc
     WRITE(*,'(7x,"                  lsf : ",i0)') lsf
     WRITE(*,'(7x,"Cutoff kinetic energy : ",e12.5)') ecutwfc
     !
  END IF
  !
  CALL mp_bcast(nq1,        ionode_id, world_comm )
  CALL mp_bcast(nq2,        ionode_id, world_comm )
  CALL mp_bcast(nq3,        ionode_id, world_comm )
  CALL mp_bcast(start_q,    ionode_id, world_comm )
  CALL mp_bcast(last_q,     ionode_id, world_comm )
  CALL mp_bcast(nci,        ionode_id, world_comm )
  CALL mp_bcast(laddxc,     ionode_id, world_comm )
  CALL mp_bcast(lsf,        ionode_id, world_comm )
  CALL mp_bcast(ecutwfc,    ionode_id, world_comm )
  !
  ! Compute irreducible q grid
  !
  nqbz = nq1 * nq2 * nq3
  fildyn       = 'matdyn'
  CALL q_points()
  IF(last_q == 0) last_q = nqs
  !
  RETURN
  !
100 WRITE(*,*) "Stop in stdin. reading namelist file"
  !
  WRITE(*,'(7x,"               q grid : ",3(i0,2x))') nq1, nq2, nq3
  WRITE(*,'(7x," # of Matsubara freq. : ",i0)') nci
  WRITE(*,'(7x,"          The first q : ",i0)') start_q
  WRITE(*,'(7x,"           The last q : ",i0)') last_q
  WRITE(*,'(7x,"               laddxc : ",i0)') laddxc
  WRITE(*,'(7x,"Cutoff kinetic energy : ",e12.5)') ecutwfc
  !
  CALL errore ('stdin_Kel', 'reading namelist Kel', 1)
  !
END SUBROUTINE stdin_Kel
!
! Standard input
!
SUBROUTINE stdin_scdft()
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_bcast
  USE io_global, ONLY : ionode, ionode_id
  USE control_flags, ONLY : niter, tr2
  USE input_parameters, ONLY : electron_maxstep, conv_thr
  USE wvfct, ONLY : nbnd
  USE output,        ONLY : fildyn
  USE constants, ONLY: K_BOLTZMANN_RY
  !
  USE sctk_val, ONLY : beta, emax, emin, fbee, lbee, ne, nmf, nx, xic, mf, wmf, &
  &                    zero_kelvin, lsf
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  REAL(dp) :: temp
  LOGICAL :: spin_fluc
  !
  NAMELIST /scdft/ temp, fbee, lbee, xic, nmf, nx, ne, emin, emax, &
  &                electron_maxstep, conv_thr, fildyn, spin_fluc
  !
  IF(ionode) THEN
     !
     temp = 0.0_dp
     fbee = 1
     lbee = nbnd
     xic = -1.0_dp
     nmf = 10
     nx = 100
     ne = 50
     emin = 1.0e-7_dp
     emax = 5.0_dp
     electron_maxstep = 100
     conv_thr = 1.0e-15_dp
     fildyn       = 'matdyn'
     spin_fluc = .FALSE.
     !
     READ(5,scdft,err=100)
     !
     WRITE(*,'(7x,"             Temparature[K] : ",e12.5)') temp
     IF (temp < 1.0e-10_dp) THEN
        zero_kelvin = .TRUE.
        beta = 1.0_dp
     ELSE
        zero_kelvin = .FALSE.
        beta = 1.0_dp / (temp*K_BOLTZMANN_RY)
        WRITE(*,'(7x,"   Inverse temparature[/Ry] : ",e12.5)') beta
     END IF
     WRITE(*,'(7x,"                  Xi cutoff : ",e12.5)') xic
     WRITE(*,'(7x,"                 First band : ",i0)') fbee
     WRITE(*,'(7x,"                  Last band : ",i0)') lbee
     WRITE(*,'(7x," # of Matsubara frequencies : ",i0)') nmf
     WRITE(*,'(7x,"                    # of xi : ",i0)') nx
     WRITE(*,'(7x,"     # of energy for QP-DOS : ",i0)') ne
     WRITE(*,'(7x,"       Minimum energy scale : ",e12.5)') emin
     WRITE(*,'(7x,"Max energy for QP-DOS [meV] : ",e12.5)') emax
     WRITE(*,'(7x,"      Convergense threshold : ",e12.5)') conv_thr
     WRITE(*,'(7x,"               Max itration : ",i0)') electron_maxstep
     WRITE(*,'(7x,"           Spin-fluvtuation : ",l)') spin_fluc
     IF(spin_fluc) THEN
        lsf = 2
     ELSE
        lsf = 1
     END IF
     !
  END IF
  !
  CALL mp_bcast(beta,             ionode_id, world_comm )
  CALL mp_bcast(zero_kelvin,      ionode_id, world_comm )
  CALL mp_bcast(xic,              ionode_id, world_comm )
  CALL mp_bcast(fbee,             ionode_id, world_comm )
  CALL mp_bcast(lbee,             ionode_id, world_comm )
  CALL mp_bcast(nmf,              ionode_id, world_comm )
  CALL mp_bcast(nx,               ionode_id, world_comm )
  CALL mp_bcast(ne,               ionode_id, world_comm )
  CALL mp_bcast(emin,             ionode_id, world_comm )
  CALL mp_bcast(emax,             ionode_id, world_comm )
  CALL mp_bcast(electron_maxstep, ionode_id, world_comm )
  CALL mp_bcast(conv_thr,         ionode_id, world_comm )
  CALL mp_bcast(lsf,              ionode_id, world_comm )
  !
  niter = electron_maxstep
  tr2 = conv_thr
  !
  ! Set Matsubara frequency and weight for the integration
  !
  ALLOCATE(mf(nmf), wmf(nmf,2))
  !
  CALL weightspoints_gl(nmf,mf,wmf(1:nmf,2))
  !
  wmf(1:nmf,1) = 2.0_dp / (pi * (mf(1:nmf)**2 + 1.0_dp)) * wmf(1:nmf,2)
  wmf(1:nmf,2) = 2.0_dp * mf(1:nmf) / (pi * (mf(1:nmf)**2 + 1.0_dp)**2) * wmf(1:nmf,2)
  mf(1:nmf) = (1.0_dp + mf(1:nmf)) / (1.0_dp - mf(1:nmf))
  !
  RETURN
  !
100 WRITE(*,*) "stop. reading namelist file"
  !
  CALL errore ('stdin_scdft', 'reading namelist SCDFT', 1)
  !
END SUBROUTINE stdin_scdft
  !
END MODULE sctk_stdin
