!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM sctk_eliashberg
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, ry_to_kelvin
  !
  IMPLICIT NONE
  !
  INTEGER :: ne, iter, maxiter, ie, nmf, imf, jmf, ii
  REAL(DP) :: temp, mix, omega_ph, lambda, mu_0, mu_c, ef, ms, omega_pl, emax, dosf, &
  &           res_chi, res_phi, res_Z, thr = 1.0e-8_dp, c1, c2, dxi, dmf, Vuc, Wint
  REAL(DP),ALLOCATABLE :: chi(:), phi(:), Z(:), mf(:), &
  &                       chi_new(:), phi_new(:), Z_new(:), &
  &                       chi0(:,:), xi(:), wght(:,:,:), dos(:), wght_esum(:,:)
  !
  Vuc = 110.5871_dp ! a.u.^3
  ef = 0.824_dp ! Ry
  ms = 1.05 !
  dosf = 2.74834984056_dp ! /Ry/cell/spin
  lambda = 0.402_dp * 2.0_dp
  mu_0 = 0.25070_dp
  mu_c = 0.556686_dp ! approx
  omega_pl = 1.5_dp ! approx
  omega_ph = 302.0_dp / ry_to_kelvin ! K -> Ry
  emax = 5.0_dp ! Ry
  !
  maxiter = 100
  mix = 0.3_dp
  temp = 5.0_dp
  temp = temp / ry_to_kelvin ! K -> Ry
  nmf = 5000
  ne = 200
  !
  ALLOCATE(mf(nmf), chi(nmf), phi(nmf), Z(nmf), &
  &        chi_new(nmf), phi_new(nmf), Z_new(nmf), &
  &        chi0(nmf,ne), xi(ne), wght(2,nmf,ne), dos(ne), wght_esum(2,nmf))
  !
  ! Matsubara frequency
  !
  DO imf = 1, nmf
    mf(imf) = REAL(2*(imf-1-nmf/2)+1, DP) * pi * temp
  END DO
  !
  ! DOS of free electron
  !
  DO ie = 1, ne
    xi(ie) = emax * REAL(ie-1,DP) / REAL(ne,DP) - ef
    dos(ie) = SQRT(ms**3*(xi(ie)+ef)) / pi**2 * 0.25_dp * Vuc
  END DO
  !
  Z= 1.0_dp
  chi = 1.0_dp
  phi = 1.0_dp
  OPEN(200, file="eliashberg.dat")
  DO iter = 1, maxiter
    !
    DO ie = 1, ne
      chi0(1:nmf,ie) = (xi(ie) + chi(1:nmf)) / SQRT((mf(1:nmf)*Z(1:nmf))**2+phi(1:nmf)**2)
    END DO
    !
    ! Weight from the trapezoidal rule
    !
    wght(1:2,1:nmf,1:ne) = 0.0_dp
    DO ie = 1, ne - 1
      DO imf = 1, nmf
        !
        c1 = chi0(imf,ie)
        c2 = chi0(imf,ie+1)
        dxi = xi(ie+1) - xi(ie)
        !
        IF(ABS(c1 - c2) < thr) THEN
          wght(1,imf,ie  ) = wght(1,imf,ie  ) + dxi * 0.5_dp / (1+c1**2)
          wght(1,imf,ie+1) = wght(1,imf,ie+1) + dxi * 0.5_dp / (1+c2**2)
          wght(2,imf,ie  ) = wght(1,imf,ie  ) + dxi * 0.5_dp * c1 / (1+c1**2)
          wght(2,imf,ie+1) = wght(1,imf,ie+1) + dxi * 0.5_dp * c2 / (1+c2**2)
        ELSE
          wght(1,imf,ie  ) = wght(1,imf,ie  ) &
          &                + dxi*(c2*(ATAN(c2)-ATAN(c1))/(c2-c1)-0.5_dp*LOG((1.0_dp+c2**2)/(1.0_dp+c1**2))/(c2-c1))/(c2-c1)
          wght(1,imf,ie+1) = wght(1,imf,ie+1) &
          &                + dxi*(c1*(ATAN(c1)-ATAN(c2))/(c1-c2)-0.5_dp*LOG((1.0_dp+c1**2)/(1.0_dp+c2**2))/(c1-c2))/(c1-c2)
          !
          wght(2,imf,ie  ) = wght(1,imf,ie  ) &
          &                + dxi*((ATAN(c2)-ATAN(c1))/(c2-c1)+0.5_dp*c2*LOG((1.0_dp+c2**2)/(1.0_dp+c1**2))/(c2-c1))/(c2-c1)
          wght(2,imf,ie+1) = wght(1,imf,ie+1) &
          &                + dxi*((ATAN(c1)-ATAN(c2))/(c1-c2)+0.5_dp*c1*LOG((1.0_dp+c1**2)/(1.0_dp+c2**2))/(c1-c2))/(c1-c2)
        END IF
        !
      END DO
    END DO
    !
    DO imf = 1, nmf
      DO ii = 1, 2
        wght_esum(ii,imf) = SUM(wght(ii,imf,1:ne) * dos(1:ne))
      END DO
    END DO
    !
    ! Compute the Left-hand side of the Eliashberg eq.
    !
    Z_new(  1:nmf) = 1.0_dp
    chi_new(1:nmf) = 0.0_dp
    phi_new(1:nmf) = 0.0_dp
    DO imf = 1, nmf
      !
      DO jmf = 1, nmf
        dmf = mf(imf) - mf(jmf)
        Wint = (omega_pl**2*mu_0+dmf**2*mu_c)/(omega_pl**2+dmf**2) - omega_ph**2*lambda/(omega_ph**2+dmf**2)
        Wint = Wint /dosf
        !
        Z_new(imf) = Z_new(imf) - temp / mf(imf) * Wint * mf(jmf) * Z(jmf) / (mf(jmf)**2*Z(jmf)**2+phi(jmf)**2)*wght_esum(1,jmf)
        chi_new(imf) = chi_new(imf) + temp / SQRT((mf(jmf)*Z(jmf))**2+phi(jmf)**2)*wght_esum(2,jmf)
        phi_new(imf) = phi_new(imf) - temp * Wint * phi(jmf) / (mf(jmf)**2*Z(jmf)**2+phi(jmf)**2)*wght_esum(1,jmf)
      END DO ! jmf
      !
    END DO ! imf
    chi_new(1:nmf) = 0.0_dp ! Ignore this level-shift term
    !
    ! Plain mixing
    !
    res_Z   = SQRT(SUM((Z_new(  1:nmf)-Z(  1:nmf))**2)) / REAL(nmf)
    res_chi = SQRT(SUM((chi_new(1:nmf)-chi(1:nmf))**2)) / REAL(nmf)
    res_phi = SQRT(SUM((phi_new(1:nmf)-phi(1:nmf))**2)) / REAL(nmf)
    !
    WRITE(*,*) " Residual Z, chi, phi", iter, res_Z, res_chi, res_phi
    !
    DO imf = 1, nmf
      WRITE(200,*) iter, mf(imf), Z(imf), phi(imf), chi(imf)
    END DO
    WRITE(200,*) ""
    !
    Z(  1:nmf) = mix*Z_new(  1:nmf) + (1.0_dp-mix)*Z(  1:nmf)
    chi(1:nmf) = mix*chi_new(1:nmf) + (1.0_dp-mix)*chi(1:nmf)
    phi(1:nmf) = mix*phi_new(1:nmf) + (1.0_dp-mix)*phi(1:nmf)
    !
  END DO ! iter
  !
  CLOSE(200)
  !
END PROGRAM sctk_eliashberg
