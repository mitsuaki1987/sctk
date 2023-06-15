!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_ini_delta
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Define initial delta
!
SUBROUTINE ini_delta(lallocate)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_bcast, mp_max
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE modes, ONLY : nmodes
  USE disp,  ONLY : nqs
  USE sctk_val, ONLY : bindx, delta, dk, dx0, kindx, emin, &
  &                     ngap, ngap1, ngap2, nx, omg0, xi, xi0
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(IN) :: lallocate
  !
  INTEGER :: it, fi = 10, is, cnt, dsp
  REAL(dp) :: thr, Z0, dosf
  CHARACTER(1) :: tmp
  !
  IF(ionode) OPEN(fi, file = 'delta.dat',status="old", action = 'read',iostat = is)
  !
  CALL mp_bcast(is,    ionode_id, world_comm )
  !
  IF(is == 0) THEN
     !
     IF(ionode) THEN
        !
        READ(fi,*) tmp, ngap1, ngap2
        !
        ngap = MAX(ngap1, ngap2)
        WRITE(*,'(7x,"Number of total points for gap equation : ",2(i0,2x))') ngap1, ngap2
        IF(lallocate) &
        &  ALLOCATE(xi(ngap,2), delta(ngap,2), dk(ngap,2), kindx(ngap,2), bindx(ngap,2))
        xi(   1:ngap,1:2) = 0.0_dp
        delta(1:ngap,1:2) = 0.0_dp
        dk(   1:ngap,1:2) = 0.0_dp
        kindx(1:ngap,1:2) = 0
        bindx(1:ngap,1:2) = 0
        !
        DO it = 1, ngap1
           READ(fi,*) xi(it,1), delta(it,1), Z0, dk(it,1), kindx(it,1), bindx(it,1)
        END DO
        !
        DO it = 1, ngap2
           READ(fi,*) xi(it,2), delta(it,2), Z0, dk(it,2), kindx(it,2), bindx(it,2)
        END DO
        !
        close(fi)
        !
     END IF ! (ionode)
     !
     CALL mp_bcast(ngap,  ionode_id, world_comm )
     CALL mp_bcast(ngap1, ionode_id, world_comm )
     CALL mp_bcast(ngap2, ionode_id, world_comm )
     IF((.NOT. ionode) .AND. lallocate) &
     &  ALLOCATE(xi(ngap,2), delta(ngap,2), dk(ngap,2), &
     &                    kindx(ngap,2), bindx(ngap,2))
     CALL mp_bcast(xi,    ionode_id, world_comm )
     CALL mp_bcast(delta, ionode_id, world_comm )
     CALL mp_bcast(dk,    ionode_id, world_comm )
     CALL mp_bcast(kindx, ionode_id, world_comm )
     CALL mp_bcast(bindx, ionode_id, world_comm )
     !
     IF(MAXVAL(ABS(delta(1:ngap1,1))) > 1e-8_dp) THEN
        !
        WRITE(stdout,'(7x,"Initial delta is read from file (delta.dat)")')
        !
        dosf = SUM(dk(1:ngap1,1), ABS(xi(1:ngap1,1)) < emin)
        dosf = dosf / dx0(minloc(ABS(xi0(1:nx)), 1))
        WRITE(stdout,'(7x,"DOS(E_F)[Ry^-1/cell/spin] : ",e12.5)') dosf
        !
        dosf = SUM(dk(1:ngap2,2), ABS(xi(1:ngap2,2)) < emin)
        dosf = dosf / dx0(minloc(ABS(xi0(1:nx)), 1))
        WRITE(stdout,'(7x,"DOS(E_F)[Ry^-1/cell/spin] : ",e12.5)') dosf
        !
        GOTO 5
        !
     END IF ! (MAXVAL(ABS(delta(1:ngap1,1))) > 1e-8_dp)
     !
     ! If delta read from file is too small
     !
     WRITE(stdout,'(7x,"Delta from file is too small !")')
     !
     DEALLOCATE(xi, delta, dk, kindx, bindx)
     !
  END IF ! (is == 0)
  !
  WRITE(stdout,'(7x,"Initial delta is theta function")')
  !
  CALL compute_d3k()
  !
  ALLOCATE(delta(ngap,2))
  !
  CALL cnt_and_dsp(nqs,cnt,dsp)
  thr = MAXVAL(omg0(1:nmodes,1:cnt))
  CALL mp_max(thr, world_comm)
  !
  IF(ionode) THEN
     CALL random_seed()
     CALL random_number(delta(1:ngap,1:2))
     delta(1:ngap,1:2) = delta(1:ngap,1:2) * thr
  END IF
  CALL mp_bcast(delta,    ionode_id, world_comm )
  !
  !delta(1:ngap,1:2) = thr**2 / (thr + xi(1:ngap,1:2))
  !
5 continue
  !
END SUBROUTINE ini_delta
!
! Compute energy groid
!
SUBROUTINE energy_grid()
  !
  USE kinds, ONLY : DP
  USE klist, ONLY : nks
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE wvfct, ONLY : et
  USE ener, ONLY : ef
  USE sctk_val, ONLY : dx0, emin, nx, xi0, Zemin
  !
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  INTEGER :: ix, imin
  REAL(dp) :: xmax, xmin, xx, ww, rhs(nx), rnx
  !
  ! Make logscale
  !
  xmax = MAXVAL(et(elph_nbnd_min:elph_nbnd_max,1:nks) - ef)
  xmin = MINVAL(et(elph_nbnd_min:elph_nbnd_max,1:nks) - ef)
  !
  ALLOCATE(xi0(nx), dx0(nx))
  CALL weightspoints_gl(nx,xi0,dx0)
  !
  rnx = REAL(nx, dp)
  !
  DO ix = 1, nx
     rhs(ix) = (1.0_dp + xi0(ix) - 2.0_dp / rnx) &
     &         * LOG(  xmax / (0.5_dp * rnx * emin * (1.0_dp - xi0(ix)))) &
     &       - (1.0_dp - xi0(ix) - 2.0_dp / rnx) &
     &         * LOG(- xmin / (0.5_dp * rnx * emin * (1.0_dp + xi0(ix))))
     rhs(ix) = ABS(rhs(ix))
  END DO
  !
  imin = MINLOC(rhs(1:nx),1)
  xx = xi0(imin)
  ww = MAX(LOG(- xmin / (0.5_dp * rnx * emin * (1.0_dp + xx))) / (1.0_dp + xx - 2.0_dp / rnx), &
  &        LOG(  xmax / (0.5_dp * rnx * emin * (1.0_dp - xx))) / (1.0_dp - xx - 2.0_dp / rnx)  )
  !
  DO ix = 1, nx
     dx0(ix) = dx0(ix) * (1.0_dp + ww * ABS(xi0(ix) - xx)) &
     &                        * emin * 0.5_dp * rnx * exp(ww * (ABS(xi0(ix) - xx) - 2.0_dp / rnx)) 
     xi0(ix) = (xi0(ix) - xx) * emin * 0.5_dp * rnx * exp(ww * (ABS(xi0(ix) - xx) - 2.0_dp / rnx))
  END DO
  emin = ABS(dx0(imin+1)) * 0.5_dp
  Zemin = ABS(dx0(imin)) * 0.5_dp
  !
END SUBROUTINE energy_grid
!
! Cp,@ite DOS_k
!
SUBROUTINE compute_dosk(dos)
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE start_k, ONLY : nk1, nk2, nk3
  USE disp,  ONLY : nq1, nq2, nq3
  USE klist, ONLY : nks
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE sctk_val, ONLY : nqbz, nx
  !
  USE sctk_tetra, ONLY : calc_dosk, interpol_indx
  IMPLICIT NONE
  !
  REAL(dp),INTENT(OUT) :: dos(nx*(elph_nbnd_max-elph_nbnd_min+1),nqbz,2)
  !
  INTEGER :: ik, ikv(3), kintp(8)
  REAL(dp) :: kv(3), wintp(1,8), &
  &            dosd(nx*(elph_nbnd_max-elph_nbnd_min+1),1,nks)
  !
  CALL calc_dosk(dosd)
  !
  ! Interpolation of weight
  !
  dos(1:nx*(elph_nbnd_max-elph_nbnd_min+1),1:nqbz,1:2) = 0.0_dp
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks,nqbz,nx,dosd,dos,elph_nbnd_min,elph_nbnd_max,nk1,nk2,nk3,nq1,nq2,nq3) &
  !$OMP & PRIVATE(ik,ikv,kv,wintp,kintp)
  !
  !$OMP DO REDUCTION(+: dos)
  DO ik = 1, nks
     !
     ikv(1) = (ik - 1) / (nk2*nk3)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
     dos(            1:nx*(elph_nbnd_max-elph_nbnd_min+1),kintp(1:8),1) &
     & = dos(        1:nx*(elph_nbnd_max-elph_nbnd_min+1),kintp(1:8),1) &
     & + MATMUL(dosd(1:nx*(elph_nbnd_max-elph_nbnd_min+1),1:1,ik), wintp(1:1,1:8))
     !
     ikv(1) = (ik - 1) / (nk2*nk3)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp) - 0.5_dp / REAL((/nq1,nq2,nq3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
     dos(            1:nx*(elph_nbnd_max-elph_nbnd_min+1),kintp(1:8),2) &
     & = dos(        1:nx*(elph_nbnd_max-elph_nbnd_min+1),kintp(1:8),2) &
     & + MATMUL(dosd(1:nx*(elph_nbnd_max-elph_nbnd_min+1),1:1,ik), wintp(1:1,1:8))
     !
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  !
  CALL mp_sum( dos, world_comm )
  !
  CALL symm_dosk(dos)
  !
END SUBROUTINE compute_dosk
!
! Compute integration weights
!
SUBROUTINE compute_d3k()
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE sctk_val, ONLY : bindx, dk, fbee, kindx, lbee, &
  &                     nqbz, ngap, ngap1, ngap2, nx, xi, xi0
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, ix
  REAL(dp) :: dos(nx,elph_nbnd_min:elph_nbnd_max,nqbz,2), thr = 1e-12_dp, shift
  !
  CALL compute_dosk(dos)
  !
  ix = minloc(ABS(xi0(1:nx)),1)
  WRITE(stdout,'(7x,"DOS(E_F)[/Ry/cell/spin] : ",e12.5)') SUM(dos(ix,elph_nbnd_min:elph_nbnd_max,1:nqbz,1)) 
  WRITE(stdout,'(7x,"DOS(E_F)[/Ry/cell/spin] : ",e12.5)') SUM(dos(ix,elph_nbnd_min:elph_nbnd_max,1:nqbz,2)) 
  !
  ! Query of ngap1
  !
  ngap1 = 0
  DO ik = 1, nqbz
     DO ib = fbee, lbee
        !
        IF(elph_nbnd_min <= ib .AND. ib <= elph_nbnd_max) THEN
           ngap1 = ngap1 + MAX(1, count(dos(1:nx,ib,ik,1) > thr))
        ELSE
           ngap1 = ngap1 + 1
        END IF
        !
     END DO
  END DO
  !
  ! Query of ngap2
  !
  ngap2 = 0
  DO ik = 1, nqbz
     DO ib = fbee, lbee
        !
        IF(elph_nbnd_min <= ib .AND. ib <= elph_nbnd_max) THEN
           ngap2 = ngap2 + MAX(1, count(dos(1:nx,ib,ik,2) > thr))
        ELSE  
           ngap2 = ngap2 + 1
        END IF
        !
     END DO
  END DO
  !
  WRITE(stdout,'(7x,"Number of total points for gap equation : ",2(i0,2x))') ngap1, ngap2
  ngap = MAX(ngap1, ngap2)
  ALLOCATE(xi(ngap,2), dk(ngap,2), kindx(ngap,2), bindx(ngap,2))
  !
  xi(   1:ngap,1:2) = 0.0_dp
  dk(   1:ngap,1:2) = 0.0_dp
  kindx(1:ngap,1:2) = 0
  bindx(1:ngap,1:2) = 0
  !
  ! Map xi, d3k, dosk, indx, bindx
  !
  shift = 0.0_dp
  CALL store_d3k(dos(1:nx,elph_nbnd_min:elph_nbnd_max,1:nqbz,1),shift,ngap1,xi(1:ngap1,1), &
  &              dk(1:ngap1,1),kindx(1:ngap1,1),bindx(1:ngap1,1))
  shift = 0.5_dp
  CALL store_d3k(dos(1:nx,elph_nbnd_min:elph_nbnd_max,1:nqbz,2),shift,ngap2,xi(1:ngap2,2), &
  &              dk(1:ngap2,2),kindx(1:ngap2,2),bindx(1:ngap2,2))
  !
END SUBROUTINE compute_d3k
!
! Map d3k, dosk, xi, kindx, bindx
!
SUBROUTINE store_d3k(dos,shift,ngap,xi,dk,kindx,bindx)
  !
  USE kinds, ONLY : DP
  USE start_k, ONLY : nk1, nk2, nk3
  USE disp,  ONLY : nq1, nq2, nq3
  USE wvfct, ONLY : et
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE ener, ONLY : ef
  !
  USE sctk_val, ONLY : dx0, fbee, lbee, nqbz, nx, xi0
  USE sctk_tetra, ONLY : interpol_indx
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ngap
  REAL(dp),INTENT(IN) :: dos(nx,elph_nbnd_min:elph_nbnd_max,nqbz), shift
  REAL(dp),INTENT(OUT) :: xi(ngap), dk(ngap)
  INTEGER,INTENT(OUT) :: kindx(ngap), bindx(ngap)
  !
  INTEGER :: ngap0, ik, ib, ix, kintp(8), ikv(3)
  REAL(dp) :: thr = 1e-12_dp, kv(3), wintp(8)
  !
  ngap0 = 0
  DO ik = 1, nqbz
     !
     ikv(1) = (ik - 1) / (nq2*nq3)
     ikv(2) = (ik - 1 - ikv(1)*nq2*nq3) / nq3
     ikv(3) =  ik - 1 - ikv(1)*nq2*nq3 - ikv(2)*nq3
     !
     DO ib = fbee, lbee
        !
        IF(elph_nbnd_min <= ib .AND. ib <= elph_nbnd_max) THEN
           !
           IF(COUNT(dos(1:nx,ib,ik) > thr) >= 2) THEN
              !
              DO ix = 1, nx                        
                 !
                 IF(dos(ix,ib,ik) > thr) THEN
                    !
                    ngap0 = ngap0 + 1
                    !
                    xi(ngap0) = xi0(ix)
                    dk(ngap0) = dos(ix,ib,ik) * dx0(ix)
                    kindx(ngap0) = ik
                    bindx(ngap0) = ib
                    !
                 END IF
                 !
              END DO ! ix
              !
              CYCLE
              !
           END IF
           !
        END IF
        !
        ! If this state is far from Fermi level
        !
        ngap0 = ngap0 + 1
        !
        kv(1:3) = (REAL(ikv(1:3), dp) + shift) / REAL((/nq1,nq2,nq3/), dp)
        CALL interpol_indx((/nk1,nk2,nk3/),kv,kintp,wintp)
        xi(ngap0) = DOT_PRODUCT(wintp(1:8), et(ib,kintp(1:8))) - ef
        !
        dk(ngap0) = 1.0_dp / REAL(nqbz, dp)
        kindx(ngap0) = ik
        bindx(ngap0) = ib
        !
     END DO ! ib
     !
  END DO ! ik
  !
  IF(ngap0 /= ngap) THEN
     CALL errore ('store_d3k', 'ngap0 /= ngap', ngap)
  END IF
  !
END SUBROUTINE store_d3k
!
! Symmetrize Dos
!
SUBROUTINE symm_dosk(dos)
  !
  USE kinds, ONLY : DP
  USE symm_base, ONLY : nsym, s
  USE disp,  ONLY : nq1, nq2, nq3
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE sctk_val, ONLY : nqbz, nx
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(inout) :: dos(nx,elph_nbnd_min:elph_nbnd_max,nqbz,2)
  !
  INTEGER :: ik, ik2, isym, ikv2(3), nksym(nqbz)
  REAL(dp) :: dostmp(nx,elph_nbnd_min:elph_nbnd_max,nqbz), kv1(3), kv2(3)
  !
  ! 1
  !
  nksym(1:nqbz) = 0
  dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,1:nqbz) = 0.0_dp
  !
  DO ik = 1, nqbz
     !
     ikv2(1) = (ik - 1) / (nq2*nq3)
     ikv2(2) = (ik - 1 - ikv2(1)*nq2*nq3) / nq3
     ikv2(3) =  ik - 1 - ikv2(1)*nq2*nq3 - ikv2(2)*nq3
     kv1(1:3) = REAL(ikv2(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
     !
     DO isym = 1, nsym
        !
        kv2(1:3) = MATMUL(REAL(s(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nq1,nq2,nq3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1e-8_dp)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nq1,nq2,nq3/))
        ik2 = 1 + ikv2(3) + ikv2(2)*nq3 + ikv2(1)*nq2*nq3
        !
        nksym(ik2) = nksym(ik2) + 1
        dostmp(    1:nx,elph_nbnd_min:elph_nbnd_max,ik2) &
        & = dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,ik2) &
        & + dos(   1:nx,elph_nbnd_min:elph_nbnd_max,ik,1)
        !
     END DO
     !
  END DO
  !
  DO ik = 1, nqbz
     dos(       1:nx,elph_nbnd_min:elph_nbnd_max,ik,1) &
     & = dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,ik) / REAL(nksym(ik), dp)
  END DO
  !
  ! 2
  !
  nksym(1:nqbz) = 0
  dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,1:nqbz) = 0.0_dp
  !
  DO ik = 1, nqbz
     !
     ikv2(1) = (ik - 1) / (nq2*nq3)
     ikv2(2) = (ik - 1 - ikv2(1)*nq2*nq3) / nq3
     ikv2(3) =  ik - 1 - ikv2(1)*nq2*nq3 - ikv2(2)*nq3
     kv1(1:3) = (REAL(ikv2(1:3), dp) + 0.5_dp) / REAL((/nq1,nq2,nq3/), dp)
     !
     DO isym = 1, nsym
        !
        kv2(1:3) = MATMUL(REAL(s(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nq1,nq2,nq3/), dp)
        kv2(1:3) = kv2(1:3) - 0.5_dp
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1e-8_dp)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nq1,nq2,nq3/))
        ik2 = 1 + ikv2(3) + ikv2(2)*nq3 + ikv2(1)*nq2*nq3
        !
        nksym(ik2) = nksym(ik2) + 1
        dostmp(    1:nx,elph_nbnd_min:elph_nbnd_max,ik2) &
        & = dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,ik2) &
        & + dos(   1:nx,elph_nbnd_min:elph_nbnd_max,ik,2)
        !
     END DO
     !
  END DO
  !
  DO ik = 1, nqbz
     dos(       1:nx,elph_nbnd_min:elph_nbnd_max,ik,2) &
     & = dostmp(1:nx,elph_nbnd_min:elph_nbnd_max,ik) / REAL(nksym(ik), dp)
  END DO
  !
END SUBROUTINE symm_dosk
!
! Initialization for lambda and mu
!
SUBROUTINE ini_lambda_mu()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE io_global, ONLY : stdout
  USE start_k, ONLY : nk1, nk2, nk3
  USE klist, ONLY : nks
  USE disp,  ONLY : nq1, nq2, nq3
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE fermisurfer_common, ONLY : b_low, b_high
  USE sctk_val, ONLY : bindx, dk, kindx, nqbz, ngap, nx, xi0, dx0
  !
  USE sctk_tetra, ONLY : calc_dosk, interpol_indx
  !
  IMPLICIT NONE
  INTEGER :: ik, ikv(3), ib, kintp(8)
  REAL(dp) :: dos(elph_nbnd_min:elph_nbnd_max,nqbz,2), &
  &          dosd(elph_nbnd_min:elph_nbnd_max,1,nks), kv(3), wintp(1,8)
  !
  nx = 1
  DEALLOCATE(xi0, dx0)
  ALLOCATE(xi0(nx), dx0(nx))
  xi0(1) = 0d0
  dx0(1) = 1d0
  !
  CALL calc_dosk(dosd)
  !
  ! Interpolation of weight
  !
  dos(elph_nbnd_min:elph_nbnd_max,1:nqbz,1:2) = 0.0_dp
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nqbz,nx,dosd,dos,elph_nbnd_min,elph_nbnd_max, &
  !$OMP &        nk1,nk2,nk3,nq1,nq2,nq3,nks) &
  !$OMP & PRIVATE(ik,ikv,kv,wintp,kintp)
  !
  !$OMP DO REDUCTION(+: dos)
  DO ik = 1, nks
     !
     ikv(1) = (ik - 1) / (nk2*nk3)
     ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
     ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
     kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp) - 0.5_dp / REAL((/nq1,nq2,nq3/), dp)
     CALL interpol_indx((/nq1,nq2,nq3/),kv,kintp,wintp)
     dos(            elph_nbnd_min:elph_nbnd_max,             kintp(1:8),2) &
     & = dos(        elph_nbnd_min:elph_nbnd_max,             kintp(1:8),2) &
     & + MATMUL(dosd(elph_nbnd_min:elph_nbnd_max,1:1,ik), wintp(1:1,1:8))
     !
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  !
  CALL mp_sum( dos, world_comm )
  !
  CALL symm_dosk(dos)
  !
  WRITE(stdout,'(7x,"DOS(E_F)[/Ry/cell/spin] : ",e12.5)') SUM(dos(elph_nbnd_min:elph_nbnd_max,1:nqbz,2)) 
  !
  ngap = nqbz * (b_high-b_low+1)
  WRITE(stdout,'(7x,"Number of total points for gap equation : ",i0)') ngap
  ALLOCATE(dk(ngap,1), kindx(ngap,1), bindx(ngap,1))
  !
  ngap = 0
  DO ik = 1, nqbz
     DO ib = b_low, b_high
        ngap = ngap + 1
        dk(ngap,1) = dos(ib,ik,2)
        kindx(ngap,1) = ik
        bindx(ngap,1) = ib
     END DO
  END DO
  !
END SUBROUTINE ini_lambda_mu
!
END MODULE sctk_ini_delta
