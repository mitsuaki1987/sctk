!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_gapeq_rhs
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Make Effective interaction
!
SUBROUTINE make_effint()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : world_comm
  USE modes, ONLY : nmodes
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : beta, bindx, effint, gg, kindx, ngap, &
  &                    omg, Vc, xi, nci, freq_min
  !
  USE sctk_kernel_weight, ONLY : Kweight, calc_Kel
  !
  IMPLICIT NONE
  !
  INTEGER :: igap, ik, ib, jgap, jk, jb, imode, ngap1(2)
  REAL(dp) :: x, xp, om, Kph(2), Kel
  !
  CALL start_clock("make_effint")
  !
  CALL divide(world_comm, ngap(1), ngap1(1), ngap1(2))
  IF(.NOT. ALLOCATED(effint)) ALLOCATE(effint(ngap(2),ngap1(1):ngap1(2),2))
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ngap1,ngap,nmodes,kindx,bindx,beta,nci, &
  !$OMP &        xi,omg,gg,Vc,effint,elph_nbnd_min,elph_nbnd_max,freq_min) &
  !$OMP & PRIVATE(igap,jgap,ik,jk,ib,jb,imode,x,xp,om,Kph,Kel)
  !
  !$OMP DO
  DO igap = ngap1(1), ngap1(2)
     !
     x = xi(igap,1)
     ik = kindx(igap,1)
     ib = bindx(igap,1)
     !
     DO jgap = 1, ngap(2)
        !
        xp = xi(jgap,2)
        jk = kindx(jgap,2)
        jb = bindx(jgap,2)
        !
        ! Electron-phonon
        !
        Kph(1:2) = 0.0_dp
        !
        IF(ALL((/elph_nbnd_min <= ib, ib <= elph_nbnd_max, &
        &        elph_nbnd_min <= jb, jb <= elph_nbnd_max/))) THEN
           !
           DO imode = 1, nmodes
              !
              om = omg(imode,jk,ik)
              IF(om < freq_min) CYCLE
              !
              Kph(1) = Kph(1) + gg(imode,jb,jk,ib,ik) * Kweight(beta, ABS(x), ABS(xp), om)
              Kph(2) = Kph(2) + gg(imode,jb,jk,ib,ik) * Kweight(beta, ABS(xp), ABS(x), om)
              !       
           END DO ! imode
           !
        END IF
        !
        ! Coulomb (+ spin-fluctuation)
        !
        Kel = calc_Kel(ABS(x)+ABS(xp),Vc(1:nci,jb,jk,ib,ik))
        !
        effint(jgap,igap,1:2) = Kph(1:2) + Kel
        !
     END DO ! jgap
     !
  END DO ! igap
  !$OMP END DO
  !$OMP END PARALLEL
  !
  CALL stop_clock("make_effint")
  !
END SUBROUTINE make_effint
!
! Make RHS vector for linear probrem.
!
SUBROUTINE gapeq_rhs(delta,res)
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE constants, ONLY : RYTOEV
  !
  USE sctk_val, ONLY : beta, zero_kelvin, dk, effint, &
  &                    ngapmax, ngap, xi, xic, Z, dxq
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(in)  :: delta(ngapmax,2)
  REAL(dp),INTENT(OUT) :: res(ngapmax,2)
  !
  INTEGER :: igap, ngap1(2), nh(2), ii
  REAL(dp) :: chi(ngapmax,2), dosh(2), Kh, dlth(2), xmax(2), tanhd, tanhe
  !
  CALL start_clock("gapeq_rhs")
  !
  CALL divide(world_comm, ngap(1), ngap1(1), ngap1(2))
  !
  res(1:ngapmax,1:2) = 0.0_dp
  chi(1:ngapmax,1:2) = 0.0_dp
  !
  IF(zero_kelvin) THEN
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & SHARED(ngap,xi,delta,dxq,chi) &
    !$OMP & PRIVATE(igap,ii)
    !
    DO ii = 1, 2
      !
      !$OMP DO
      DO igap = 1, ngap(ii)
        IF(xi(igap,ii)**2 + delta(igap,ii)**2 > dxq(igap,ii)**2) THEN
          chi(igap,ii) = 1.0_dp
        ELSE
          chi(igap,ii) = 0.0_dp
        END IF
      END DO ! igap
      !$OMP END DO NOWAIT
      !
    END DO ! ii
    !
    !$OMP END PARALLEL
  ELSE
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & SHARED(ngap,xi,delta,dxq,chi,beta) &
    !$OMP & PRIVATE(igap,tanhd,tanhe,ii)
    !
    DO ii = 1, 2
      !
      !$OMP DO
      DO igap = 1, ngap(ii)
        !
        tanhe = TANH(0.5_dp * beta * SQRT(xi(igap,ii)**2 + delta(igap,ii)**2))
        tanhd = TANH(0.5_dp * beta * ABS(dxq(igap,ii)))
        !
        IF(1.0_dp - tanhd * tanhe < 1.0e-10_dp) THEN
          IF(xi(igap,ii)**2 + delta(igap,ii)**2 > dxq(igap,ii)**2) THEN
            chi(igap,ii) = tanhe * (1.0_dp + tanhd) / (1.0_dp + tanhd * tanhe)
          ELSE
            chi(igap,ii) = 0.0_dp
          END IF
        ELSE
          chi(igap,ii) = tanhe * (1.0_dp + tanhd) / (1.0_dp + tanhd * tanhe) &
          &                    * (1.0_dp - tanhd) / (1.0_dp - tanhd * tanhe)
        END IF
        !
      END DO ! igap
      !$OMP END DO NOWAIT
      !
    END DO ! ii
    !
    !$OMP END PARALLEL
  END IF
  !
  DO ii = 1, 2
    chi(1:ngap(ii),ii) = 0.5_dp * dk(1:ngap(ii),ii) * delta(1:ngap(ii),ii) * chi(1:ngap(ii),ii) &
    &                  / SQRT(xi(1:ngap(ii),ii)**2 + delta(1:ngap(ii),ii)**2)
  END DO
  !
  CALL dgemv("T", ngap(2), ngap1(2)-ngap1(1)+1, 1.0_dp, effint(1:ngap(2),ngap1(1):ngap1(2),1), ngap(2), &
  &    chi(1:ngap(2),2), 1, 1.0_dp, res(ngap1(1):ngap1(2),1), 1)
  !
  CALL dgemv("N", ngap(2), ngap1(2)-ngap1(1)+1, 1.0_dp, effint(1:ngap(2),ngap1(1):ngap1(2),2), ngap(2), &
  &    chi(ngap1(1):ngap1(2),1), 1, 1.0_dp, res(1:ngap(2),2), 1)
  !
  ! High energy region
  !
  IF(xic > 0) THEN
    !
    DO ii = 1, 2
      nh(ii) = count(xi(1:ngap(ii),ii) > xic)
      dlth(ii) = SUM(delta(1:ngap(ii),ii) / xi(1:ngap(ii),ii), xi(1:ngap(ii),ii) > xic) &
      &        / SUM(1.0_dp / xi(1:ngap(ii),ii)**2,         xi(1:ngap(ii),ii) > xic)
      !
      xmax(ii) = MAXVAL(xi(1:ngap(ii),ii))
      dosh(ii) = SUM(dk(1:ngap(ii),ii), xi(1:ngap(ii),ii) > xic) / (xmax(ii) - xic)
      !
      WRITE(stdout,'(9x,"Extrapol ",i1," : ",e12.5)') ii, dlth(ii) / xic * RYTOEV * 1.0e3_dp
    END DO
    !
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & SHARED(res,ngap,effint,ngap1,xi,xic,xmax,nh,dosh,dlth) &
    !$OMP & PRIVATE(igap,Kh)
    !$OMP DO
    DO igap = 1, ngap(2)
      Kh = SUM(effint(igap,ngap1(1):ngap1(2),2), xi(ngap1(1):ngap1(2),1) > xic) / REAL(nh(1), dp)
      res(igap,2) = res(igap,2) + 0.5_dp * dosh(1) * Kh * dlth(1) / xmax(1)
    END DO
    !$OMP END DO
    !     
    !$OMP DO
    DO igap = ngap1(1), ngap1(2)
      Kh = SUM(effint(1:ngap(2),igap,1), xi(1:ngap(2),2) > xic) / REAL(nh(2), dp)
      res(igap,1) = res(igap,1) + 0.5_dp * dosh(2) * Kh * dlth(2) / xmax(2)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    !
  END IF ! (xic > 0)
  !
  CALL mp_sum( res, world_comm )
  !
  res(1:ngapmax,1:2) = - res(1:ngapmax,1:2) / (1.0_dp + Z(1:ngapmax,1:2)) - delta(1:ngapmax,1:2)
  !
  CALL stop_clock("gapeq_rhs")
  !
END SUBROUTINE gapeq_rhs
!
! Make Kel, Kph & a part of Jacobian
!
SUBROUTINE gapeq_rhs_qpdos()
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE fermisurfer_common,   ONLY : b_low, b_high
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE sctk_val, ONLY : beta, zero_kelvin, bindx, delta, dk, dltF, ggf, kindx, &
  &                    nci, ngap, nx, omgf, VcF, xi, xi0, xic, ZF, freq_min
  !
  USE sctk_kernel_weight, ONLY : Kweight, calc_Kel
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, imode, ix, nh, nks0, nks1
  REAL(dp) :: x, xp, om, eqp, Kph, Kel, Kave, dlth, dosh, xmax
  !
  CALL start_clock("gapeq_rhs_qpdos")
  !
  ALLOCATE(dltF(nx,b_low:b_high,nks))
  dltF(1:nx,b_low:b_high,1:nks) = 0.0_dp
  !
  ! Extrapolation of high energy region
  !
  nh = count(xi(1:ngap(2),2) > xic)
  dlth = SUM(delta(1:ngap(2),2) / xi(1:ngap(2),2), xi(1:ngap(2),2) > xic) &
  &    / SUM(1.0_dp / xi(1:ngap(2),2)**2,      xi(1:ngap(2),2) > xic)
  !
  xmax = MAXVAL(xi(1:ngap(2),2))
  dosh = SUM(dk(1:ngap(2),2), xi(1:ngap(2),2) > xic) / (xmax - xic)
  !
  WRITE(stdout,*) nh, dlth, xmax, dosh
  !
  CALL divide(world_comm, nks,nks0, nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(xi0,nx,nks0,nks1,ngap,b_low,b_high,nmodes,omgf,ggf,VcF, &
  !$OMP &        xi,dk,delta,kindx,bindx,dltF,ZF,beta,elph_nbnd_min,elph_nbnd_max, &
  !$OMP &        nh,dlth,dosh,xic,xmax,zero_kelvin,nci,freq_min) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,imode,ix, &
  !$OMP &         x,xp,om,eqp,Kph,Kel,Kave)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        DO ix = 1, nx 
           !
           Kave = 0.0_dp
           x = ABS(xi0(ix))
           !
           DO jgap = 1, ngap(2)
              !
              xp = ABS(xi(jgap,2))
              jk = kindx(jgap,2)
              jb = bindx(jgap,2)
              !
              ! Electron-phonon
              !
              Kph = 0.0_dp
              !
              IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
                 !
                 DO imode = 1, nmodes
                    !
                    om = omgf(imode,jk,ik)
                    IF(om < freq_min) CYCLE
                    !
                    Kph = Kph + ggf(imode,jb,jk,ib,ik) * Kweight(beta, ABS(x), ABS(xp), om)
                    !
                 END DO ! imode
                 !               
              END IF ! (elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max)
              !
              ! Coulomb (+ spin-fluctuation)
              !
              Kel = calc_Kel(ABS(x)+ABS(xp),VcF(1:nci,jb,jk,ib,ik))
              !
              ! Add RHS of the gap eq.
              !
              eqp = SQRT(delta(jgap,2)**2 + xp**2)
              !
              IF(zero_kelvin) THEN
                 dltF(ix,ib,ik) = dltF(ix,ib,ik) &
                 &              - 0.5_dp * dk(jgap,2) * (Kph + Kel) / (1.0_dp + ZF(ix,ib,ik)) &
                 &              * delta(jgap,2) / eqp
              ELSE
                 dltF(ix,ib,ik) = dltF(ix,ib,ik) &
                 &              - 0.5_dp * dk(jgap,2) * (Kph + Kel) / (1.0_dp + ZF(ix,ib,ik)) &
                 &              * delta(jgap,2) / eqp * TANH(0.5_dp * beta * eqp)
              END IF
              !
              IF(xi(jgap,2) > xic) Kave = Kave + (Kph + Kel)
              !
           END DO ! jgap
           !
           ! High-energy extrapolation
           !
           IF(xic > 0) THEN
              !
              Kave = Kave / REAL(nh, dp)
              dltF(ix,ib,ik) = dltF(ix,ib,ik) - 0.5_dp * Kave / (1.0_dp + ZF(ix,ib,ik)) &
              &                          * dosh * dlth / xmax
              !
           END IF ! (xic > 0)
           !
           ! Delta -> quasi-particle energy
           !
           dltF(ix,ib,ik) = SQRT(xi0(ix)**2 + dltF(ix,ib,ik)**2)
           !
        END DO ! ix
        !
     END DO ! ib
     !
  END DO ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  ! 
  CALL mp_sum( dltF, world_comm )
  !
  DEALLOCATE(ggf, VcF, omgf)
  !
  CALL stop_clock("gapeq_rhs_qpdos")
  !
END SUBROUTINE gapeq_rhs_qpdos
!
! Make Kel, Kph & a part of Jacobian
!
SUBROUTINE gapeq_rhs_f()
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE fermisurfer_common,   ONLY : b_low, b_high
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE constants, ONLY : RYTOEV
  !
  USE sctk_val, ONLY : beta, bindx, delta, dk, dltf, ggf, kindx, freq_min, &
  &                    nci, ngap, omgf, VcF, xi, xic, Zf, zero_kelvin
  !
  USE sctk_kernel_weight, ONLY : Kweight, calc_Kel
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, imode, nh, nks0, nks1
  REAL(dp) :: xp, om, eqp, Kph, Kel, Kave, dlth, dosh, xmax
  !
  CALL start_clock("gapeq_rhs_f")
  !
  ALLOCATE(dltf(1,b_low:b_high,nks))
  dltf(1,b_low:b_high,1:nks) = 0.0_dp
  !
  ! Extrapolation of high energy region
  !
  nh = count(xi(1:ngap(2),2) > xic)
  dlth = SUM(delta(1:ngap(2),2) / xi(1:ngap(2),2), xi(1:ngap(2),2) > xic) &
  &    / SUM(1.0_dp / xi(1:ngap(2),2)**2,        xi(1:ngap(2),2) > xic)
  !
  xmax = MAXVAL(xi(1:ngap(2),2))
  dosh = SUM(dk(1:ngap(2),2), xi(1:ngap(2),2) > xic) / (xmax - xic)
  !
  IF(xic > 0.0_dp) WRITE(stdout,*) dlth / xic * RYTOEV * 1.0e3_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,ngap,b_low,b_high,nmodes,omgf,ggf,VcF,nci, &
  !$OMP &        xi,dk,delta,kindx,bindx,dltf,Zf,beta,elph_nbnd_min,elph_nbnd_max, &
  !$OMP &        nh,dlth,dosh,xic,xmax,zero_kelvin,freq_min) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,imode,xp,om,eqp,Kph,Kel,Kave)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        Kave = 0.0_dp
        !
        DO jgap = 1, ngap(2)
           !
           xp = ABS(xi(jgap,2))
           jk = kindx(jgap,2)
           jb = bindx(jgap,2)
           !
           ! Electron-phonon
           !
           Kph = 0.0_dp
           !
           IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
              !
              DO imode = 1, nmodes
                 !
                 om = omgf(imode,jk,ik)
                 IF(om < freq_min) CYCLE
                 !
                 Kph = Kph + ggf(imode,jb,jk,ib,ik) * Kweight(beta, 0.0_dp, ABS(xp), om)
                 !
              END DO ! imode
              !
           END IF
           !
           ! Coulomb (+ spin-fluctuation)
           !
           Kel = calc_Kel(ABS(xp),VcF(1:nci,jb,jk,ib,ik))
           !
           ! Add RHS of the gap eq.
           !
           eqp = SQRT(delta(jgap,2)**2 + xp**2)
           !
           IF(zero_kelvin) THEN
              dltf(1,ib,ik) = dltf(1,ib,ik) &
              &             - 0.5_dp * dk(jgap,2) * (Kph + Kel) / (1.0_dp + Zf(1,ib,ik)) &
              &             * delta(jgap,2) / eqp
           ELSE
              dltf(1,ib,ik) = dltf(1,ib,ik) &
              &             - 0.5_dp * dk(jgap,2) * (Kph + Kel) / (1.0_dp + Zf(1,ib,ik)) &
              &             * delta(jgap,2) / eqp * TANH(0.5_dp * beta * eqp)
           END IF
           !
           IF(xi(jgap,2) > xic) Kave = Kave + (Kph + Kel)
           !
        END DO ! jgap
        !
        ! High-energy extrapolation
        !
        IF(xic > 0.0_dp) THEN
           !
           Kave = Kave / REAL(nh, dp)
           !
           dltf(1,ib,ik) = dltf(1,ib,ik) - 0.5_dp * Kave / (1.0_dp + Zf(1,ib,ik)) &
           &                          * dosh * dlth / xmax
           !
        END IF ! (xic > 0)
        !
     END DO ! ib
     !
  END DO ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  CALL mp_sum( dltf, world_comm )
  !
  ! Ry -> meV
  !
  dltf(1,b_low:b_high,1:nks) = dltf(1,b_low:b_high,1:nks) * 13605.692283_dp
  !
  DEALLOCATE(ggf, VcF, omgf)
  !
  CALL stop_clock("gapeq_rhs_f")
  !
END SUBROUTINE gapeq_rhs_f
!
! Calc lambda_k and mu_k
!
SUBROUTINE make_lambda_mu_f()
  !
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE sctk_kernel_weight, ONLY : calc_Kel
  USE kinds, ONLY : DP
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE fermisurfer_common,   ONLY : b_low, b_high
  USE io_global, ONLY : stdout
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : bindx, dk, ggf, kindx, ngapmax, omgf, ZF, dltf, VcF, nci, freq_min
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  USE sctk_tetra, ONLY : calc_dosk
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, nks0, nks1
  REAL(dp) :: lambda, mu, weight, dosd(elph_nbnd_min:elph_nbnd_max,nks)
  !
  CALL start_clock("make_lambda_mu_f")
  !
  ALLOCATE(ZF(1,b_low:b_high,nks), dltf(1,b_low:b_high,nks)) !dltf -> mu_k
  ZF(1,b_low:b_high,1:nks) = 0.0_dp
  dltf(1,b_low:b_high,1:nks) = 0.0_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,b_low,b_high,ngapmax,kindx,bindx,ZF,dltf,dk,ggf,omgf,VcF,nmodes,nci,freq_min) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        DO jgap = 1, ngapmax
           !
           jk = kindx(jgap,1)
           jb = bindx(jgap,1)
           !
           ZF(1,ib,ik) = ZF(1,ib,ik) &
           &         + 2.0_dp * dk(jgap,1) * SUM(ggf(1:nmodes,jb,jk,ib,ik) / ABS(omgf(1:nmodes,jk,ik)), &
           &                                     omgf(1:nmodes,jk,ik) > freq_min)
           !
           dltf(1,ib,ik) = dltf(1,ib,ik) + dk(jgap,1) * calc_Kel(0.0_dp,VcF(1:nci,jb,jk,ib,ik))
           !
        END DO ! jgap
        !
     END DO ! ib
     !
  END DO ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  CALL calc_dosk(dosd)
  !
  CALL mp_sum( ZF, world_comm )
  CALL mp_sum( dltf, world_comm )
  CALL mp_sum( dosd, world_comm )
  !
  lambda = SUM(ZF(1,b_low:b_high,1:nks)*dosd(b_low:b_high,1:nks))
  mu = SUM(  dltf(1,b_low:b_high,1:nks)*dosd(b_low:b_high,1:nks))
  weight = SUM(                         dosd(b_low:b_high,1:nks))
  lambda = lambda / weight
  mu = mu / weight
  !
  WRITE(stdout,'(7x,"    mu : ",e15.5)') mu
  WRITE(stdout,'(7x,"lambda : ",e15.5)') lambda
  !
  CALL stop_clock("make_lambda_mu_f")
  !
END SUBROUTINE make_lambda_mu_f
!
END MODULE sctk_gapeq_rhs
