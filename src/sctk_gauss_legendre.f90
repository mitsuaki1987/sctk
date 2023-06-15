!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_gauss_legendre
  !
  IMPLICIT NONE
  !
  public weightspoints_gl
  private legendre
  !
CONTAINS
!
! Weights & Points for Gauss-Legendre method
!
SUBROUTINE weightspoints_gl(n,x,w)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY: pi
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  REAL(dp),INTENT(OUT) :: x(n), w(n)
  !
  INTEGER :: i, itr
  REAL(dp) :: Pn, Pnm1, PP
  !
  Pn = 0.0_dp
  Pnm1 = 0.0_dp
  !
  DO i = 1, n
     !
     x(i) = COS(pi * (REAL(n - i + 1, dp) - 0.25_dp) / (REAL(n, dp) + 0.5_dp))
     !
     DO itr = 1, 100
        !
        CALL legendre(n    , x(i), Pn  )
        CALL legendre(n - 1, x(i), Pnm1)
        !
        PP = REAL(n, dp) * (Pnm1 - x(i) * Pn) / (1.0_dp - x(i)**2)
        !
        x(i) = x(i) - Pn / PP
        !
        IF(ABS(Pn/PP) < 1e-10_dp) EXIT
        !
     END DO
     !
     IF(itr >= 100) THEN
        CALL errore ('weightspoints_gl', 'newton.', n)
     END IF
     !
     CALL legendre(n - 1, x(i), Pnm1)     
     w(i) = 2.0_dp * (1.0_dp - x(i)**2) / (REAL(n, dp) * Pnm1)**2
     !
  END DO
  !
END SUBROUTINE weightspoints_gl
!
! Legendre polynomial
!
SUBROUTINE legendre(n,x,P)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  REAL(dp),INTENT(IN) :: x
  REAL(dp),INTENT(OUT) :: P
  !
  INTEGER :: i
  REAL(dp) :: Pold, Pnew
  !
  IF(n == 0) THEN
     P = 1.0_dp
     RETURN
  ELSE IF(n == 1) THEN
     P = x
     RETURN
  END IF
  !
  Pold = 1.0_dp
  Pnew = x
  !
  DO i = 1, n - 1
     !
     P = (REAL(2 * i + 1, dp) * x * Pnew - REAL(i, dp) * Pold) / REAL(i + 1, dp)
     Pold = Pnew
     Pnew = P
     !
  END DO
  !
END SUBROUTINE legendre
!
END MODULE sctk_gauss_legendre
