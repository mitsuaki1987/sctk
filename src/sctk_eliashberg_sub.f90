!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_eliashberg_sub
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP),SAVE :: temp, omega_pl, omega_ph, lambda, mu_0, mu_c, dosf
  REAL(DP),ALLOCATABLE,SAVE :: dos(:), xi(:)
  !
  CONTAINS
  !
SUBROUTINE eliashberg_rhs(Zcp,rZcp)
  !
  USE kinds, ONLY : DP
  USE sctk_val, ONLY : nmf, ne, mf
  !
  REAL(DP),INTENT(IN) :: Zcp(nmf,3) ! Z, chi, phi
  REAL(DP),INTENT(OUT) :: rZcp(nmf,3) ! OUT - IN
  !
  INTEGER :: ie, imf, jmf, ii
  REAL(DP) :: chi0(nmf,ne), wght(2,nmf,ne), c1, c2, dxi, wght_esum(2,nmf), dmf, &
  &           thr = 1.0e-8_dp, Wint
  !
  DO ie = 1, ne
    chi0(1:nmf,ie) = (xi(ie) + Zcp(1:nmf,2)) / SQRT((mf(1:nmf)*Zcp(1:nmf,1))**2+Zcp(1:nmf,3)**2)
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
  rZcp(1:nmf,1) = 1.0_dp ! Z
  rZcp(1:nmf,2) = 0.0_dp ! chi 
  rZcp(1:nmf,3) = 0.0_dp ! phi
  DO imf = 1, nmf
    !
    DO jmf = 1, nmf
      dmf = mf(imf) - mf(jmf)
      Wint = (omega_pl**2*mu_0+dmf**2*mu_c)/(omega_pl**2+dmf**2) - omega_ph**2*lambda/(omega_ph**2+dmf**2)
      Wint = Wint /dosf
      !
      rZcp(imf,1) = rZcp(imf,1) - temp / mf(imf) * Wint * mf(jmf) * Zcp(jmf,1) &
      &                         / (mf(jmf)**2*Zcp(jmf,1)**2+Zcp(jmf,3)**2)*wght_esum(1,jmf)
      rZcp(imf,2) = rZcp(imf,2) + temp / SQRT((mf(jmf)*Zcp(jmf,1))**2+Zcp(jmf,3)**2)*wght_esum(2,jmf)
      rZcp(imf,3) = rZcp(imf,3) - temp * Wint * Zcp(jmf,3) / (mf(jmf)**2*Zcp(jmf,1)**2+Zcp(jmf,3)**2)*wght_esum(1,jmf)
    END DO ! jmf
    !
  END DO ! imf
  rZcp(1:nmf,2) = 0.0_dp ! Ignore this level-shift term
  !
  rZcp(1:nmf,1:3) = rZcp(1:nmf,1:3) - Zcp(1:nmf,1:3)
  !
END SUBROUTINE eliashberg_rhs
  !
END MODULE sctk_eliashberg_sub