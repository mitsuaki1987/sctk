!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!> @mainpage Density Functional theory for superconductors
!>
!> This package contains the following programs.
!>
!> ### sctk.x ( scdft_main.f90 ) ###
!>
!> Compute @f$\Delta_{n k}@f$ by solving the Kohn-Sham gap equation.
!> Compute electron-electron Coulomb interaction via RPA
!> It perform a post process 
!> Compute @f$\Delta_{n k}@f$ on a dense @f$k@f$ grid with non-scf
!> calculation and  output them for the fermisurfer
!> It perform a post process 
!> Compute the quasi particle DOS 
!> @f[
!>  D(\varepsilon) = \sum_{n k} \delta \left(\varepsilon_{n k} - 
!>  \sqrt{\xi_{n k}^2 + |\Delta_{n k}|^2} \right)
!> @f]
!> It perform a post process of scdft.x
!> Compute the ultrasonic attenuation coefficient 
!> 
!
!> SCDFT\@QE Main routine to solve the gap equation
!> @f[
!>   \Delta_{n k} = - \frac{1}{2} \sum_{n' k'} K_{n k n' k'} 
!>   \frac{\Delta_{n' k'}}{E_{n' k'}}
!>   \tanh \left( \frac{\beta E_{n' k'}}{2} \right)
!> @f]
!> and obtain @f$\Delta_{n k}@f$ .
!>
!> @author Mitsuaki Kawamura
!> 
PROGRAM sctk_main
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE environment,ONLY : environment_start, environment_end
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE mp_pools, ONLY : npool
  USE input_parameters, ONLY : calculation, restart_mode
  USE mp_world, ONLY : mpime, nproc
  USE mp_exx, ONLY : negrp
  USE wvfct, ONLY : et
  USE control_ph, ONLY : start_q, last_q, current_iq
  USE elph_tetra_mod, ONLY : lshift_q
  USE constants, ONLY: eV_to_kelvin, K_BOLTZMANN_RY
  USE exx, ONLY : exx_fft_create
  USE uspp, ONLY : okvan
  USE control_flags,  ONLY : tqr
  USE io_files, ONLY : prefix, tmp_dir
  !
  USE sctk_wfc, ONLY : get_wfcg, fft_wfc
  USE sctk_dmuxc, ONLY : generate_dmuxc, apply_xc
  USE sctk_coulomb, ONLY : make_scrn, make_kel, Kel_frequency, Coulomb_parameter, &
  &                       prepare_q, chebyshev_interpol, write_Kel
  USE sctk_invert, ONLY : invert
  USE sctk_io_delta, ONLY : read_delta, write_dos, output_frmsf
  USE sctk_stdin, ONLY : stdin_scdft, stdin_Kel, stdin_control
  USE sctk_read_file, ONLY : read_elph, read_Coulomb, read_a2Fsave, degenerated_band
  USE sctk_rotate_kernel, ONLY : expand_g_v, expand_g_v_f
  USE sctk_ini_delta, ONLY : ini_delta, ini_lambda_mu, energy_grid
  USE sctk_z, ONLY : make_Z, make_Z_f, make_Z_qpdos
  USE sctk_gapeq_rhs, ONLY : make_effint, gapeq_rhs_f, &
  &                          gapeq_rhs_qpdos, make_lambda_mu_f
  USE sctk_broyden, ONLY : broyden_gapeq, compute_dabs
  !USE wvfct, ONLY : nbnd ! debug
  !USE sctk_val, ONLY : Kel, nk_p, k0_p, nbnd_p, bnd0_p, nmf, nqbz ! debug
  !USE mp, ONLY : mp_sum ! debug
  !USE mp_world, ONLY : world_comm ! debug
  USE sctk_val, ONLY : dltF, ZF, laddxc, Wscr, lsf, beta, zero_kelvin, &
  &                    bisec_min, bisec_max, bisec_step, delta, ngap
  USE sctk_usonic, ONLY : calc_fvel, calc_usonic
  USE sctk_qpdos, ONLY : egrid, calc_sdos
  USE sctk_spinfluc, ONLY : lambda_sf
  USE sctk_clock, ONLY : sctk_print_clock
  USE sctk_io_delta, ONLY : out_delta
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, ibisec
  REAL(dp) :: dabs, delta0, tcmid, tcmin, tcmax
  LOGICAL :: needwf = .FALSE.
  !REAL(DP),ALLOCATABLE :: Kel_all(:,:,:,:,:) !debug
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'SCTK' )
  CALL sctk_write_hash()
  !
  IF (nproc /= npool*negrp) &
  & CALL errore ('sctk_main', 'npool(-nk) times band-group(-nb) must equal to the number of processes.', npool)
  !
  CALL stdin_control()
  CALL read_file_new ( needwf )
  caLL degenerated_band()
  CALL read_a2Fsave()
  lshift_q = .TRUE.
  !
  IF(TRIM(calculation) == "kel") THEN
     !
     WRITE(stdout,'(/,5x,"#####  Read from STDIN  #####",/)') 
     CALL stdin_Kel()
     IF(okvan) tqr = .TRUE.
     CALL exx_fft_create()
     !
     WRITE(stdout,'(/,5x,"#####  Read gvectors.dat & evc.dat  #####",/)')
     CALL get_wfcg()
     !
     WRITE(stdout,'(/,5x,"#####  FFT wfc(G) -> wfc(r)  #####",/)')
     CALL fft_wfc()
     !
     IF(laddxc>0 .OR. lsf>0) THEN
        !
        WRITE(stdout,'(/,5x,"#####  Generate f_xc  #####",/)')
        CALL generate_dmuxc()
        !
     END IF
     !
     WRITE(stdout,'(/,5x,"#####  Set frequency grid  #####",/)')
     CALL Kel_frequency()
     !
     WRITE(stdout,'(/,5x,"#####  Compute K_el  #####",/)')
     !
     DO iq = start_q, last_q
        !
        current_iq = iq
        !
        WRITE(stdout,'(/,7x,"q-point : ",i0,/)') current_iq
        CALL prepare_q()
        CALL make_scrn()
        !
        ! Coulomb part
        !
        CALL apply_xc()
        !
        IF(lsf>0) CALL lambda_sf()
        !
        CALL make_kel()
        DEALLOCATE(Wscr)
        !
        CALL Coulomb_parameter()
        !
        !!!!!!!!  start debug
        !ALLOCATE(Kel_all(0:nmf+1,nbnd,nbnd,nqbz,lsf+1))
        !Kel_all(0:nmf+1,1:nbnd,1:nbnd,1:nqbz,1:lsf+1) = 0.0_dp
        !Kel_all(0:nmf+1,1:nbnd,bnd0_p+1:bnd0_p+nbnd_p,k0_p+1:k0_p+nk_p,1:lsf+1) &
        !& = Kel(0:nmf+1,1:nbnd,       1:       nbnd_p,     1:     nk_p,1:lsf+1)
        !CALL mp_sum(Kel_all,world_comm)
        !IF(mpime==0) WRITE(80+iq,'(3e25.15)') Kel_all(:,:,:,:,1) !debug
        !IF(lsf>0) THEN
        !   IF(mpime==0) WRITE(80+iq,*) 
        !   IF(mpime==0) WRITE(80+iq,'(2e25.15)') Kel_all(0:1,:,:,:,2) !debug
        !END IF
        !DEALLOCATE(Kel_all)
        !!!!!!!!  finish debug
        CALL chebyshev_interpol()
        CALL write_Kel()
        !
     END DO ! iq = start_q, last_q
     !
  ELSE ! calculation = scdft, lambda_mu_k, qpdos, deltaf, ultrasonic
     !
     WRITE(stdout,'(/,5x,"#####  Read from STDIN  #####",/)')
     CALL stdin_scdft()
     !
     WRITE(stdout,'(/,5x,"#####  Read from elph*.dat  #####",/)')
     CALL read_elph()
     !
     WRITE(stdout,'(/,5x,"#####  Read from vc*.dat  #####",/)')
     CALL read_Coulomb()
     !
     CALL energy_grid()
     !
     IF(TRIM(calculation) == "scdft" .OR. TRIM(calculation) == "scdft_tc") THEN
        !
        WRITE(stdout,'(/,5x,"#####  Set or read initial delta  #####",/)')
        CALL ini_delta(.TRUE.)
        !
        WRITE(stdout,'(/,5x,"#####  Average matrix in grid  #####",/)')
        CALL expand_g_v()
        !
        IF(TRIM(calculation) == "scdft") THEN
           !
           WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
           CALL make_Z()
           !
           WRITE(stdout,'(/,5x,"#####  Compute effective interaction  #####",/)')
           CALL make_effint()
           !
           WRITE(stdout,'(/,5x,"#####  Solve gap equation  #####",/)')
           CALL broyden_gapeq(ngap*2,delta)
           CALL out_delta(TRIM(tmp_dir) // TRIM(prefix) // ".scgap")
           !
        ELSE
           !
           IF(bisec_min < 0.0_dp .OR. bisec_max > 0.0_dp) THEN
              !
              ! Zero kelvin
              !
              zero_kelvin = .TRUE.
              !
              WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
              CALL make_Z()
              !
              WRITE(stdout,'(/,5x,"#####  Compute effective interaction  #####",/)')
              CALL make_effint()
              !
              WRITE(stdout,'(/,5x,"#####  Solve gap equation  #####",/)')
              CALL broyden_gapeq(ngap*2,delta)
              CALL compute_dabs(.FALSE., delta0)
              CALL out_delta(TRIM(tmp_dir) // TRIM(prefix) // ".scgap")
              !
              WRITE(stdout,'(/,5x,"T[K]-D[meV] ",2e15.5/)') 0.0_dp, delta0
              IF(delta0 < 1.0e-3_dp) THEN
                 CALL errore ('sctk_main', 'No SC even at the zero Kelvin.', 1)
              ELSE
                 tcmin = 0.0_dp
                 restart_mode = "restart"
              END IF
              !
              ! Find upper limit of Tc
              !
              DO ibisec = 1, 10
                 !
                 ! BCS Tc estimation
                 ! 2 Delta / (k_B Tc) = 3.54
                 !
                 tcmax = 2.0_dp / 3.54_dp * delta0 * 1.0e-3_dp * eV_to_kelvin * REAL(ibisec, dp)
                 beta = 1.0_dp / (tcmax*K_BOLTZMANN_RY)
                 zero_kelvin = .FALSE.
                 !
                 CALL ini_delta(.FALSE.)
                 !
                 WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
                 CALL make_Z()
                 !
                 WRITE(stdout,'(/,5x,"#####  Compute effective interaction  #####",/)')
                 CALL make_effint()
                 !
                 WRITE(stdout,'(/,5x,"#####  Solve gap equation  #####",/)')
                 CALL broyden_gapeq(ngap*2,delta)
                 CALL compute_dabs(.FALSE., dabs)
                 !
                 !
                 WRITE(stdout,'(/,5x,"T[K]-D[meV]: ",2e15.5,/)') tcmax, dabs
                 !
                 IF(dabs < delta0*0.001_dp) EXIT
                 !
              END DO ! ibisec = 1, 10
              !
           END IF
           !
           ! Find Tc by the bisection method
           !
           DO ibisec = 1, bisec_step
              !
              tcmid = 0.5_dp * (tcmin + tcmax)
              beta = 1.0_dp / (tcmid*K_BOLTZMANN_RY)
              !
              CALL ini_delta(.FALSE.)
              !
              WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
              CALL make_Z()
              !
              WRITE(stdout,'(/,5x,"#####  Compute effective interaction  #####",/)')
              CALL make_effint()
              !
              WRITE(stdout,'(/,5x,"#####  Solve gap equation  #####",/)')
              CALL broyden_gapeq(ngap*2,delta)
              CALL compute_dabs(.FALSE., dabs)
              !
              WRITE(stdout,'(/,5x,"T[K]-D[meV]: ",2e15.5,/)') tcmid, dabs
              !
              IF(dabs < delta0*0.001_dp) THEN
                 tcmax = tcmid
              ELSE
                 tcmin = tcmid
              END IF
              !
           END DO ! ibisec = 1, 10
           !
           tcmid = 0.5_dp * (tcmin + tcmax)
           WRITE(stdout,'(/,5x,"Transition temperature [K] ",e15.5,/)') tcmid
           !
        END IF
        !
     ELSE ! calculation = lambda_mu_k, qpdos, deltaf, ultrasonic
        !
        WRITE(stdout,'(/,5x,"#####  Average matrix in grid  #####",/)')
        CALL expand_g_v_f()
        !
        IF(TRIM(calculation) == "lambda_mu_k") THEN
           !
           CALL ini_lambda_mu()
           !
           WRITE(stdout,'(/,5x,"#####  Integrate gapeq  #####",/)')
           CALL make_lambda_mu_f()
           !
           WRITE(stdout,'(/,5x,"#####  Write mu.frmsf and lambda_frmsf  #####",/)')
           IF(mpime == 0) THEN
              CALL output_frmsf(et, dltF, TRIM(tmp_dir) // TRIM(prefix) // "_mu.frmsf")
              CALL output_frmsf(et, ZF, TRIM(tmp_dir) // TRIM(prefix) // "_lambda.frmsf")
           END IF
           !
        ELSE ! calculation = qpdos, deltaf, ultrasonic
           !
           WRITE(stdout,'(/,5x,"#####  Read gap function  #####",/)')
           CALL read_delta()
           !
           IF(TRIM(calculation) == "qpdos") THEN
              !
              CALL egrid()
              !
              WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
              CALL make_Z_qpdos()
              !
              WRITE(stdout,'(/,5x,"#####  Integrate gapeq  #####",/)')
              CALL gapeq_rhs_qpdos()
              !
              WRITE(stdout,'(/,5x,"#####  Compute qpdos  #####",/)')
              !
              CALL calc_sdos()
              IF(mpime == 0) CALL write_dos()
              !
           ELSE ! deltaf, ultrasonic
              !
              WRITE(stdout,'(/,5x,"#####  Compute renormalization factor : Z  #####",/)')
              CALL make_Z_f()
              !
              WRITE(stdout,'(/,5x,"#####  Integrate gapeq  #####",/)')
              CALL gapeq_rhs_f()
              !
              IF(TRIM(calculation) == "deltaf") THEN
                 !
                 WRITE(stdout,'(/,5x,"#####  Write delta.frmsf  #####",/)')
                 !
                 IF(mpime == 0) THEN
                    CALL output_frmsf(et, dltF, TRIM(tmp_dir) // TRIM(prefix) // "_delta.frmsf")
                    CALL output_frmsf(et, ZF, TRIM(tmp_dir) // TRIM(prefix) // "_Z.frmsf")
                 END IF
                 !
              ELSE ! calculation = ultrasonic
                 !
                 WRITE(stdout,'(/,5x,"#####  Compute ultrasonic attenuation  #####",/)')
                 !
                 CALL calc_fvel()
                 CALL calc_usonic()
                 !
              END IF
           END IF
        END IF
     END IF
  END IF
  !
  CALL sctk_print_clock()
  !
  CALL environment_end ( 'SCTK' )
#if defined(__MPI)
  CALL mp_global_end ( )
#endif
  !
END PROGRAM
