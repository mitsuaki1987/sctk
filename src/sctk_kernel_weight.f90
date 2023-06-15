!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_kernel_weight
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Weight for Z
!
FUNCTION Zweight(x,y,z,tx,ty,tz) RESULT(Wz)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: x, y, z, tx, ty, tz
  REAL(dp) :: Wz
  !
  REAL(dp) :: thr = 1e-8_dp, txmz, txpz, mm, mp, pm, pp
  !
  IF(ABS(x) < thr) THEN
     !
     pp = 1.0_dp / (  y + z)
     mp = 1.0_dp / (- y + z)
     !
     IF(ABS(y + z) < thr) THEN
        !
        Wz = 0.125_dp * ( &
        &      (2.0_dp - 3.0_dp * ty**2) * (1.0_dp - ty**2) / (3.0_dp * ty) &
        &    + ((tz**2 - 1.0_dp) + (-2.0_dp + (tz + ty) * tz + 2.0_dp * (tz - ty) * mp) * mp / tz) * mp &
        &              )
        !
     ELSE IF(ABS(- y + z) < thr) THEN
        !
        Wz = 0.125_dp * ( &
        &      ((tz**2 - 1.0_dp) + (-2.0_dp + (tz - ty) * tz + 2.0_dp * (tz + ty) * pp) * pp / tz) * pp &
        &    + (2.0_dp - 3.0_dp * ty**2) * (1.0_dp - ty**2) / (- 3.0_dp * ty) &
        &              )
        !
     ELSE
        !
        Wz = 0.125_dp * ( &
        &      ((tz**2 - 1.0_dp) + (-2.0_dp + (tz - ty) * tz + 2.0_dp * (tz + ty) * pp) * pp / tz) * pp &
        &    + ((tz**2 - 1.0_dp) + (-2.0_dp + (tz + ty) * tz + 2.0_dp * (tz - ty) * mp) * mp / tz) * mp &
        &              )
        !
     END IF
     !
  ELSE
     !
     mm = 1.0_dp / (x - y - z)
     mp = 1.0_dp / (x - y + z)
     pm = 1.0_dp / (x + y - z)
     pp = 1.0_dp / (x + y + z)
     !
     txpz = tanh(x + z)
     txmz = tanh(x - z)
     !
     IF(ABS(x - y - z) < thr) THEN
        !
        Wz = 0.0625_dp * ( &
        &   (1.0_dp + ty) * (tx - 1.0_dp) * (1.0_dp + tz) * ty / (tz * tx) &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     ELSE IF(ABS(x - y + z) < thr) THEN
        !
        Wz = 0.0625_dp * ( &
        &   (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1.0_dp + ty) * (tx - 1.0_dp) * (1.0_dp - tz) * ty / (- tz * tx) &
        & + (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     ELSE IF(ABS(x + y - z) < thr) THEN
        !
        Wz = 0.0625_dp * ( &
        &   (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1.0_dp - ty) * (tx - 1.0_dp) * (1.0_dp + tz) * (- ty) / (tz * tx) &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     ELSE IF(ABS(x + y + z) < thr) THEN
        !
        Wz = 0.0625_dp * ( &
        &   (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1.0_dp - ty) * (tx - 1.0_dp) * (1.0_dp - tz) * (- ty) / (- tz * tx) &
        & )
        !
     ELSE 
        !
        Wz = 0.0625_dp * ( &
        &   (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1.0_dp - 1.0_dp / (  tz * tx)) * (-1.0_dp + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1.0_dp - 1.0_dp / (- tz * tx)) * (-1.0_dp + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     END IF
     !
  END IF
  !
END FUNCTION Zweight
!
! Weight for K
!
FUNCTION Kweight(x,y,z,tx,ty,tz) RESULT(Wk)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: x, y, z, tx, ty, tz
  REAL(dp) :: Wk
  !
  REAL(dp) :: thr = 1e-1_dp, zp, zm
  !
  IF(ABS(x) < thr) THEN
     !
     IF(ABS(y) < thr) THEN
        !
        Wk = ((-1.0_dp / z + 1.0_dp / tz) / z - 0.5_dp) / z
        !
     ELSE
        !
        zp = 1.0_dp / (  y + z)
        zm = 1.0_dp / (- y + z)
        !
        IF (ABS(y + z) < thr) THEN
           !
           Wk = ( &
           &        (- 1.0_dp + ty * tz + (- ty + tz) * zm) * (- zm) &
           &    ) / (4.0_dp * ty * tz)
           !
        ELSE IF (ABS(- y + z) < thr) THEN
           !
           Wk = ( &
           &        (- 1.0_dp - ty * tz + (  ty + tz) * zp) * (  zp) &
           &    ) / (4.0_dp * ty * tz)
           !
        ELSE
           !
           Wk = ( &
           &        (-1.0_dp - ty * tz + (  ty + tz) * zp) * (  zp) &
           &      + (-1.0_dp + ty * tz + (- ty + tz) * zm) * (- zm) &
           &    ) / (4.0_dp * ty * tz)
           !
        END IF
        !
     END IF
     !
  ELSE IF(ABS(y) < thr) THEN
     !
     zp = 1.0_dp / (  x + z)
     zm = 1.0_dp / (- x + z)
     !
     IF(ABS(x + z) < thr) THEN
        !
        Wk = ( &
        &      (-1.0_dp + tx * tz + (- tx + tz) * zm) * (- zm) &
        &    ) / (4.0_dp * tx * tz)
        !
     ELSE IF(ABS(- x + z) < thr) THEN
        !
        Wk = ( &
        &        (-1.0_dp - tx * tz + (  tx + tz) * zp) * (  zp) &
        &    ) / (4.0_dp * tx * tz)
        !
     ELSE
        !
        Wk = ( &
        &        (-1.0_dp - tx * tz + (  tx + tz) * zp) * (  zp) &
        &      + (-1.0_dp + tx * tz + (- tx + tz) * zm) * (- zm) &
        &    ) / (4.0_dp * tx * tz)
        !
     END IF
     !
  ELSE
     !
     IF(ABS(x - y - z) < thr) THEN
        !
        Wk =  ( &
        &         (  (1.0_dp - tx) * (  tz + 1.0_dp) * (1.0_dp + ty)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8.0_dp * tx * ty * tz)
        !
     ELSE IF(ABS(x - y + z) < thr) THEN
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- (1.0_dp - tx) * (- tz + 1.0_dp) * (1.0_dp + ty)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8.0_dp * tx * ty * tz)
        !
     ELSE IF(ABS(x + y - z) < thr) THEN
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- (1.0_dp - tx) * (  tz + 1.0_dp) * (1.0_dp - ty)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8.0_dp * tx * ty * tz)
        !
     ELSE IF(ABS(x + y + z) < thr) THEN
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  (1.0_dp - tx) * (- tz + 1.0_dp) * (1.0_dp - ty)) &
        &     ) / (8.0_dp * tx * ty * tz)
        !
     ELSE
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8.0_dp * tx * ty * tz)
        !
     END IF
     !
  END IF
  !
END FUNCTION Kweight
!
! Compute Kel
!
FUNCTION calc_Kel(x,Vc0) RESULT(Kel)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  !
  USE sctk_val, ONLY : mf, nmf, wmf, nci
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: x, Vc0(nci)
  REAL(dp) :: Kel
  !
  INTEGER :: ici
  REAL(dp) :: Vc1(nmf), mf1(nmf), x0, cheb(nmf,nci)
  !
  mf1(1:nmf) = x * mf(1:nmf)
  !
  x0 = COS(pi / REAL(2 * nci, dp))
  mf1(1:nmf) = x0 * (mf1(1:nmf) - 1.0_dp) / (mf1(1:nmf) + 1.0_dp)
  mf1(1:nmf) = ACOS(mf1(1:nmf))
  !
  DO ici = 1, nci
     cheb(1:nmf,ici) = COS(REAL(ici-1,dp) * mf1(1:nmf))
  END DO
  !
  Vc1(1:nmf) = MATMUL(cheb(1:nmf,1:nci), Vc0(1:nci))
  !
  Kel = DOT_PRODUCT(wmf(1:nmf,1), Vc1(1:nmf))
  !
END FUNCTION calc_Kel
!
! Compute Kel
!
FUNCTION calc_Zsf(x,Vc0) RESULT(Zsf)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  !
  USE sctk_val, ONLY : mf, nmf, wmf, emin, nci
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: x, Vc0(nci)
  REAL(dp) :: Zsf
  !
  INTEGER :: ici
  REAL(dp) :: Vc1(nmf), mf1(nmf), x0, cheb(nmf,nci), e0
  !
  IF(x < emin) THEN
     e0 = emin
  ELSE
     e0 = x
  END IF
  !
  mf1(1:nmf) = e0 * mf(1:nmf)
  !
  x0 = COS(pi / REAL(2 * nci, dp))
  mf1(1:nmf) = x0 * (mf1(1:nmf) - 1.0_dp) / (mf1(1:nmf) + 1.0_dp)
  mf1(1:nmf) = ACOS(mf1(1:nmf))
  !
  DO ici = 1, nci
     cheb(1:nmf,ici) = COS(REAL(ici-1,dp) * mf1(1:nmf))
  END DO
  !
  Vc1(1:nmf) = MATMUL(cheb(1:nmf,1:nci), Vc0(1:nci))
  Zsf = DOT_PRODUCT(wmf(1:nmf,2), Vc1(1:nmf)) / e0
  !
END FUNCTION calc_Zsf
!
! Weight for Z
!
FUNCTION Zweight_f(y,z,ty,tz) RESULT(Wz)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: y, z, ty, tz
  REAL(dp) :: Wz
  !
  REAL(dp) :: thr = 1e-8_dp, mp, pp
  !
  pp = 1.0_dp / (  y + z)
  mp = 1.0_dp / (- y + z)
  !
  IF(ABS(y + z) < thr) THEN
     !
     Wz = 0.125_dp * ( &
     &      (2.0_dp - 3.0_dp * ty**2) * (1.0_dp - ty**2) / (3.0_dp * ty) &
     &    + ((tz**2 - 1.0_dp) + (-2.0_dp + (tz + ty) * tz + 2.0_dp * (tz - ty) * mp) * mp / tz) * mp &
     &              )
     !
  ELSE IF(ABS(- y + z) < thr) THEN
     !
     Wz = 0.125_dp * ( &
     &      ((tz**2 - 1.0_dp) + (-2.0_dp + (tz - ty) * tz + 2.0_dp * (tz + ty) * pp) * pp / tz) * pp &
     &    + (2.0_dp - 3.0_dp * ty**2) * (1.0_dp - ty**2) / (- 3.0_dp * ty) &
     &              )
     !
  ELSE
     !
     Wz = 0.125_dp * ( &
     &      ((tz**2 - 1.0_dp) + (-2.0_dp + (tz - ty) * tz + 2.0_dp * (tz + ty) * pp) * pp / tz) * pp &
     &    + ((tz**2 - 1.0_dp) + (-2.0_dp + (tz + ty) * tz + 2.0_dp * (tz - ty) * mp) * mp / tz) * mp &
     &              )
     !
  END IF
  !
END FUNCTION Zweight_f
!
! Weight for K
!
FUNCTION Kweight_f(y,z,ty,tz) RESULT(Wk)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: y, z, ty, tz
  REAL(dp) :: Wk
  !
  REAL(dp) :: thr = 1e-4_dp, zp, zm
  !
  IF(ABS(y) < thr) THEN
     !
     Wk = ((-1.0_dp / z + 1.0_dp / tz) / z - 0.5_dp) / z
     !
  ELSE
     !
     zp = 1.0_dp / (  y + z)
     zm = 1.0_dp / (- y + z)
     !
     IF (ABS(y + z) < thr) THEN
        !
        Wk = ( &
        &        (- 1.0_dp + ty * tz + (- ty + tz) * zm) * (- zm) &
        &    ) / (4.0_dp * ty * tz)
        !
     ELSE IF (ABS(- y + z) < thr) THEN
        !
        Wk = ( &
        &        (- 1.0_dp - ty * tz + (  ty + tz) * zp) * (  zp) &
        &    ) / (4.0_dp * ty * tz)
        !
     ELSE
        !
        Wk = ( &
        &        (-1.0_dp - ty * tz + (  ty + tz) * zp) * (  zp) &
        &      + (-1.0_dp + ty * tz + (- ty + tz) * zm) * (- zm) &
        &    ) / (4.0_dp * ty * tz)
        !
     END IF
     !
  END IF
  !
END FUNCTION Kweight_f
!
END MODULE sctk_kernel_weight
