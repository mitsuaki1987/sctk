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
  USE sctk_val, ONLY : beta, zero_kelvin, bindx, effint, gg, kindx, ngap1, ngap2, omg, Vc, xi, nci
  !
  USE sctk_kernel_weight, ONLY : Kweight, calc_Kel
  !
  IMPLICIT NONE
  !
  INTEGER :: igap, ik, ib, jgap, jk, jb, im, ngap10, ngap11
  REAL(dp) :: x, xp, om, Kph, Kel, bx, bxp, tx, txp, tom
  !
  CALL start_clock("make_effint")
  !
  CALL divide(world_comm, ngap1, ngap10, ngap11)
  IF(.NOT. ALLOCATED(effint)) ALLOCATE(effint(ngap2,ngap10:ngap11))
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ngap10,ngap11,ngap1,ngap2,nmodes,kindx,bindx,beta,nci, &
  !$OMP &        xi,omg,gg,Vc,effint,elph_nbnd_min,elph_nbnd_max,zero_kelvin) &
  !$OMP & PRIVATE(igap,jgap,ik,jk,ib,jb,im,x,xp,om,Kph,Kel, &
  !$OMP &         bx,bxp,tx,txp,tom)
  !
  !$OMP DO
  DO igap = ngap10, ngap11
     !
     x = xi(igap,1)
     bx = 0.5_dp * beta * x
     tx = TANH(bx)
     ik = kindx(igap,1)
     ib = bindx(igap,1)
     !
     DO jgap = 1, ngap2
        !
        xp = xi(jgap,2)
        bxp = 0.5_dp * beta * xp
        txp = TANH(bxp)
        jk = kindx(jgap,2)
        jb = bindx(jgap,2)
        !
        ! Electron-phonon
        !
        Kph = 0.0_dp
        !
        IF(ALL((/elph_nbnd_min <= ib, ib <= elph_nbnd_max, &
        &        elph_nbnd_min <= jb, jb <= elph_nbnd_max/))) THEN
           !
           IF(zero_kelvin) THEN
              Kph = - SUM(gg(1:nmodes,jb,jk,ib,ik) / (ABS(x)+ABS(xp)+omg(1:nmodes,jk,ik)))
           ELSE
              DO im = 1, nmodes
                 !
                 om = omg(im,jk,ik) * beta * 0.5_dp
                 tom = TANH(om)
                 !
                 Kph = Kph + gg(im,jb,jk,ib,ik) * beta &
                 &         * Kweight(bx, bxp, om, tx, txp, tom)
                 !       
              END DO ! im
           END IF
           !
           Kph = Kph * 2.0_dp
           !
        END IF
        !
        ! Coulomb (+ spin-fluctuation)
        !
        Kel = calc_Kel(ABS(x)+ABS(xp),Vc(1:nci,jb,jk,ib,ik))
        !
        effint(jgap,igap) = Kph + Kel
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
SUBROUTINE gapeq_rhs(veco)
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE constants, ONLY : RYTOEV
  !
  USE sctk_val, ONLY : beta, zero_kelvin, delta, dk, effint, &
  &                    ngap, ngap1, ngap2, xi, xic, Z
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(OUT) :: veco(ngap,2)
  !
  INTEGER :: igap, ngap10, ngap11, nh
  REAL(dp) :: chi(ngap,2), dosh, Kh, dlth, xmax
  !
  CALL start_clock("gapeq_rhs")
  !
  CALL divide(world_comm, ngap1, ngap10, ngap11)
  !
  veco(1:ngap,1:2) = 0.0_dp
  chi(1:ngap,1:2) = 0.0_dp
  !
  IF(zero_kelvin) THEN
     chi(1:ngap1,1) = 1.0_dp
     chi(1:ngap2,2) = 1.0_dp
  ELSE
     chi(1:ngap1,1) = TANH(0.5_dp * beta * SQRT(xi(1:ngap1,1)**2 + delta(1:ngap1,1)**2))
     chi(1:ngap2,2) = TANH(0.5_dp * beta * SQRT(xi(1:ngap2,2)**2 + delta(1:ngap2,2)**2))
  END IF
  !
  chi(1:ngap1,1) = 0.5_dp * dk(1:ngap1,1) * delta(1:ngap1,1) * chi(1:ngap1,1) &
  &              / SQRT(xi(1:ngap1,1)**2 + delta(1:ngap1,1)**2)
  !
  chi(1:ngap2,2) = 0.5_dp * dk(1:ngap2,2) * delta(1:ngap2,2) * chi(1:ngap2,2) &
  &              / SQRT(xi(1:ngap2,2)**2 + delta(1:ngap2,2)**2)
  !
  CALL dgemv("T", ngap2, ngap11-ngap10+1, 1.0_dp, effint(1:ngap2,ngap10:ngap11), ngap2, &
  &    chi(1:ngap2,2), 1, 1.0_dp, veco(ngap10:ngap11,1), 1)
  !
  CALL dgemv("N", ngap2, ngap11-ngap10+1, 1.0_dp, effint(1:ngap2,ngap10:ngap11), ngap2, &
  &    chi(ngap10:ngap11,1), 1, 1.0_dp, veco(1:ngap2,2), 1)
  !
  ! High energy region
  !
  IF(xic > 0) THEN
     !
     ! 1, 
     !
     nh = count(xi(1:ngap1,1) > xic)
     dlth = SUM(delta(1:ngap1,1) / xi(1:ngap1,1), xi(1:ngap1,1) > xic) &
     &    / SUM(1.0_dp / xi(1:ngap1,1)**2,         xi(1:ngap1,1) > xic)
     !
     xmax = MAXVAL(xi(1:ngap1,1))
     dosh = SUM(dk(1:ngap1,1), xi(1:ngap1,1) > xic) / (xmax - xic)
     !
     WRITE(stdout,'(9x,"Extrapol 1 : ",e12.5)') dlth / xic * RYTOEV * 1.0e3_dp
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(veco,ngap2,effint,ngap10,ngap11,xi,xic,xmax,nh,dosh,dlth) &
     !$OMP & PRIVATE(igap,Kh)
     !$OMP DO
     DO igap = 1, ngap2
        Kh = SUM(effint(igap,ngap10:ngap11), xi(ngap10:ngap11,1) > xic) / REAL(nh, dp)
        veco(igap,2) = veco(igap,2) + 0.5_dp * dosh * Kh * dlth / xmax
     END DO
     !$OMP END DO
     !$OMP END PARALLEL
     !
     ! 2, 
     !
     nh = count(xi(1:ngap2,2) > xic)
     dlth = SUM(delta(1:ngap2,2) / xi(1:ngap2,2), xi(1:ngap2,2) > xic) &
     &    / SUM(1.0_dp / xi(1:ngap2,2)**2,         xi(1:ngap2,2) > xic)
     !
     xmax = MAXVAL(xi(1:ngap2,2))
     dosh = SUM(dk(1:ngap2,2), xi(1:ngap2,2) > xic) / (xmax - xic)
     !
     WRITE(stdout,'(9x,"Extrapol 2 : ",e12.5)') dlth / xic * RYTOEV * 1.0e3_dp
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(veco,ngap2,effint,ngap10,ngap11,xi,xic,xmax,nh,dosh,dlth) &
     !$OMP & PRIVATE(igap,Kh)
     !$OMP DO
     DO igap = ngap10, ngap11
        Kh = SUM(effint(1:ngap2,igap), xi(1:ngap2,2) > xic) / REAL(nh, dp)
        veco(igap,2) = veco(igap,2) + 0.5_dp * dosh * Kh * dlth / xmax
     END DO
     !$OMP END DO
     !$OMP END PARALLEL
     !
  END IF ! (xic > 0)
  !
  CALL mp_sum( veco, world_comm )
  !
  veco(1:ngap,1:2) = delta(1:ngap,1:2) + veco(1:ngap,1:2) / (1.0_dp + Z(1:ngap,1:2)) 
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
  &                    nci, ngap2, nx, omgf, VcF, xi, xi0, xic, ZF
  !
  USE sctk_kernel_weight, ONLY : Kweight, calc_Kel
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, im, ix, nh, nks0, nks1
  REAL(dp) :: x, xp, om, bx, bxp, tx, txp, tom, eqp, Kph, Kel, Kave, dlth, dosh, xmax
  !
  CALL start_clock("gapeq_rhs_qpdos")
  !
  ALLOCATE(dltF(nx,b_low:b_high,nks))
  dltF(1:nx,b_low:b_high,1:nks) = 0.0_dp
  !
  ! Extrapolation of high energy region
  !
  nh = count(xi(1:ngap2,2) > xic)
  dlth = SUM(delta(1:ngap2,2) / xi(1:ngap2,2), xi(1:ngap2,2) > xic) &
  &    / SUM(1.0_dp / xi(1:ngap2,2)**2,      xi(1:ngap2,2) > xic)
  !
  xmax = MAXVAL(xi(1:ngap2,2))
  dosh = SUM(dk(1:ngap2,2), xi(1:ngap2,2) > xic) / (xmax - xic)
  !
  WRITE(stdout,*) nh, dlth, xmax, dosh
  !
  CALL divide(world_comm, nks,nks0, nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(xi0,nx,nks0,nks1,ngap2,b_low,b_high,nmodes,omgf,ggf,VcF, &
  !$OMP &        xi,dk,delta,kindx,bindx,dltF,ZF,beta,elph_nbnd_min,elph_nbnd_max, &
  !$OMP &        nh,dlth,dosh,xic,xmax,zero_kelvin,nci) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,im,ix, &
  !$OMP &         x,xp,om,bx,bxp,tx,txp,tom,eqp,Kph,Kel,Kave)
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
           bx = ABS(0.5_dp * beta * x)
           tx = TANH(bx)
           !
           DO jgap = 1, ngap2
              !
              xp = ABS(xi(jgap,2))
              bxp = ABS(0.5_dp * beta * xp)
              txp = TANH(bxp)
              jk = kindx(jgap,2)
              jb = bindx(jgap,2)
              !
              ! Electron-phonon
              !
              Kph = 0.0_dp
              !
              IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
                 !
                 IF(zero_kelvin) THEN
                    Kph = - SUM(ggf(1:nmodes,jb,jk,ib,ik) / (ABS(x)+ABS(xp)+omgf(1:nmodes,jk,ik)))
                 ELSE
                    DO im = 1, nmodes
                       !
                       om = ABS(omgf(im,jk,ik) * 0.5_dp * beta)
                       tom = TANH(om)
                       !
                       Kph = Kph + ggf(im,jb,jk,ib,ik) &
                       &   * beta * Kweight(bx, bxp, om, tx, txp, tom)
                       !
                    END DO ! im
                 END IF
                 !               
                 Kph = Kph * 2.0_dp
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
  USE sctk_val, ONLY : beta, bindx, delta, dk, dltf, ggf, kindx, &
  &                    nci, ngap2, omgf, VcF, xi, xic, Zf, zero_kelvin
  !
  USE sctk_kernel_weight, ONLY : Kweight_f, calc_Kel
  USE sctk_gauss_legendre, ONLY : weightspoints_gl
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, im, nh, nks0, nks1
  REAL(dp) :: xp, om, eqp, Kph, Kel, Kave, dlth, dosh, xmax, bxp, txp, tom
  !
  CALL start_clock("gapeq_rhs_f")
  !
  ALLOCATE(dltf(1,b_low:b_high,nks))
  dltf(1,b_low:b_high,1:nks) = 0.0_dp
  !
  ! Extrapolation of high energy region
  !
  nh = count(xi(1:ngap2,2) > xic)
  dlth = SUM(delta(1:ngap2,2) / xi(1:ngap2,2), xi(1:ngap2,2) > xic) &
  &    / SUM(1.0_dp / xi(1:ngap2,2)**2,        xi(1:ngap2,2) > xic)
  !
  xmax = MAXVAL(xi(1:ngap2,2))
  dosh = SUM(dk(1:ngap2,2), xi(1:ngap2,2) > xic) / (xmax - xic)
  !
  IF(xic > 0.0_dp) WRITE(stdout,*) dlth / xic * RYTOEV * 1.0e3_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,ngap2,b_low,b_high,nmodes,omgf,ggf,VcF,nci, &
  !$OMP &        xi,dk,delta,kindx,bindx,dltf,Zf,beta,elph_nbnd_min,elph_nbnd_max, &
  !$OMP &        nh,dlth,dosh,xic,xmax,zero_kelvin) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,im,bxp,txp,tom, &
  !$OMP &         xp,om,eqp,Kph,Kel,Kave)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        Kave = 0.0_dp
        !
        DO jgap = 1, ngap2
           !
           xp = ABS(xi(jgap,2))
           jk = kindx(jgap,2)
           jb = bindx(jgap,2)
           bxp = 0.5_dp * beta * xp
           txp = TANH(bxp)
           !
           ! Electron-phonon
           !
           Kph = 0.0_dp
           !
           IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
              !
              IF(zero_kelvin) THEN
                 Kph = - SUM(ggf(1:nmodes,jb,jk,ib,ik) / (ABS(xp)+omgf(1:nmodes,jk,ik)))
              ELSE
                 DO im = 1, nmodes
                    !
                    om = ABS(0.5_dp * beta * omgf(im,jk,ik))
                    tom = TANH(om)
                    !
                    Kph = Kph + ggf(im,jb,jk,ib,ik) &
                    &   * beta * Kweight_f(bxp, om, txp, tom)
                    !
                 END DO ! im
              END IF
              Kph = Kph * 2.0_dp
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
  USE sctk_val, ONLY : bindx, dk, ggf, kindx, ngap, omgf, ZF, dltf, VcF, nci
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
  !$OMP & SHARED(nks0,nks1,b_low,b_high,ngap,kindx,bindx,ZF,dltf,dk,ggf,omgf,VcF,nmodes,nci) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        DO jgap = 1, ngap
           !
           jk = kindx(jgap,1)
           jb = bindx(jgap,1)
           !
           ZF(1,ib,ik) = ZF(1,ib,ik) &
           &         + 2.0_dp * dk(jgap,1) * SUM(ggf(1:nmodes,jb,jk,ib,ik) / ABS(omgf(1:nmodes,jk,ik)))
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
