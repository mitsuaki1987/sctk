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
  USE control_flags, ONLY : niter, tr2
  USE sctk_val, ONLY : nmf, ne, mf
  USE sctk_eliashberg_sub, ONLY : xi, dos, temp, omega_pl, omega_ph, lambda, mu_0, mu_c, dosf
  USE sctk_broyden, ONLY: broyden_gapeq
  USE klist, ONLY : nelec
  USE ener, ONLY : ef
  !
  IMPLICIT NONE
  !
  INTEGER :: ie, imf
  REAL(DP) :: ms, emax, Vuc, nelec1, nelec2
  REAL(DP),ALLOCATABLE :: Zcp(:,:)
  !
  Vuc = 110.5871_dp ! a.u.^3
  ef = 0.824_dp ! Ry
  ms = 1.05 !
  nelec = 3.0_dp
  dosf = 2.74834984056_dp ! /Ry/cell/spin
  lambda = 0.402_dp * 2.0_dp
  mu_0 = 0.25070_dp
  mu_c = 0.556686_dp ! approx
  omega_pl = 1.5_dp ! approx
  omega_ph = 302.0_dp / ry_to_kelvin ! K -> Ry
  emax = 5.0_dp ! Ry
  !
  niter = 100
  tr2 = 1.0e-10_dp
  temp = 5.0_dp
  temp = temp / ry_to_kelvin ! K -> Ry
  nmf = 5000
  ne = 200
  !
  ALLOCATE(mf(nmf), Zcp(nmf,3), xi(ne), dos(ne))
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
  nelec1 = 0.0_dp
  nelec2 = 0.0_dp
  DO ie = 1, ne - 1
    !
    nelec1 = nelec1 + 0.5_dp *(xi(ie+1)-xi(ie)) &
    &               * ( dos(ie  )*(1.0_dp-TANH(0.5_dp*xi(ie  )/temp)) &
    &                 + dos(ie+1)*(1.0_dp-TANH(0.5_dp*xi(ie+1)/temp)) )
    !
    nelec2 = nelec2 + 0.5_dp * (xi(ie+1)-xi(ie)) * (dos(ie  ) + dos(ie+1))
    !
    DO imf = 1, nmf
      nelec2 = nelec2 - 0.5_dp *(xi(ie+1)-xi(ie)) * 2.0_dp * temp &
      &               * ( dos(ie  )*(xi(ie  )/(mf(imf)**2+xi(ie  )**2)) &
      &                 + dos(ie+1)*(xi(ie+1)/(mf(imf)**2+xi(ie+1)**2)) )
    END DO
  END DO
  WRITE(*,*) "debug1", nelec1, nelec2
  !
  Zcp(1:nmf,1) = 1.0_dp
  Zcp(1:nmf,2) = 0.0_dp
  Zcp(1:nmf,3) = 1.0_dp
  !
  OPEN(200, file="eliashberg.dat")
  !
  CALL broyden_gapeq(3*nmf, Zcp)
  !
  CLOSE(200)
  !
END PROGRAM sctk_eliashberg
