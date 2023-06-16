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
FUNCTION Zweight(beta,x0,y0,z0) RESULT(Wz)
  !
  USE kinds, ONLY : DP
  USE sctk_val, ONLY : scdft_kernel, zero_kelvin
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: beta, x0, y0, z0
  REAL(dp) :: Wz
  !
  REAL(dp) :: thr = 1e-3_dp, tzpx, tzmx, mm, mp, pm, pp, &
  &           tx, ty, tz, x, y, z, bx, by, bz, txmy
  !
  Wz = 0.0_dp
  !
  IF(zero_kelvin) THEN
     IF(scdft_kernel == 1 .OR. scdft_kernel == 2) THEN
        Wz = -1d0/(x0+y0+z0)**2
     END IF
     !
     RETURN
     !
  END IF
  !
  bx = beta * x0
  by = beta * y0
  bz = beta * z0
  !
  IF(bx < thr) THEN
     IF(ABS(by-bz) < thr) THEN
        x = 0.0_dp
        y = bz
        z = bz
     ELSE
        x = 0.0_dp
        y = by
        z = bz
     END IF
  ELSE IF(by < thr) THEN
     IF(ABS(bx-bz) < thr) THEN
        x = bz
        y = 0.0_dp
        z = bz
     ELSE
        x = bx
        y = 0.0_dp
        z = bz
     END IF
  ELSE
     IF(ABS(bx+by-bz) < thr) THEN
        x = bx
        y = by
        z = bx + by
     ELSE IF(ABS(-bx+by+bz) < thr) THEN
        x = by + bz
        y = by
        z = bz
     ELSE IF(ABS(bx-by+bz) < thr) THEN
        x = bx
        y = bx + bz
        z = bz
     ELSE
        x = bx
        y = by
        z = bz        
     END IF
  END IF
  !
  tx = TANH(0.5_dp*x)
  ty = TANH(0.5_dp*y)
  tz = TANH(0.5_dp*z)
  mp = 1.0_dp/(-x + y + z)
  pm = 1.0_dp/( x - y + z)
  mm = 1.0_dp/(-x - y + z)
  pp = 1.0_dp/( x + y + z)
  !
  IF(scdft_kernel == 1) THEN
     !
     tzpx = TANH(0.5_dp*(z+x))
     tzmx = TANH(0.5_dp*(z-x))
     !
     IF(bx < thr) THEN
        IF(ABS(by-bz) < thr) THEN
           Wz = 0.125_dp*(((tz**2-1.0_dp)+2.0_dp*(2.0_dp*tz/z-1.0_dp)/(z*tz))/z &
           &              -(2.0_dp-3.0_dp*tz**2)*(1.0_dp-tz**2)/(3.0_dp*tz))
        ELSE
           Wz = ( &
           &      (0.25_dp*(tz**2-1.0_dp)+(-1.0_dp+(tz-ty)*tz*0.5_dp+2.0_dp*(tz+ty)*pp)*pp/tz)*pp &
           &     +(0.25_dp*(tz**2-1.0_dp)+(-1.0_dp+(tz+ty)*tz*0.5_dp+2.0_dp*(tz-ty)*pm)*pm/tz)*pm &
           &    )
        END IF
     ELSE IF(by < thr) THEN
        IF(ABS(bx-bz) < thr) THEN
           Wz = -(tx*tz + 1.0_dp)*(0.5_dp*(tzpx**2 - 1.0_dp) + tzpx*pp)*pp/(2.0_dp*tx*tz)
        ELSE
           Wz = ( &
           &      -(tx*tz + 1.0_dp)*(0.5_dp*(tzpx**2 - 1.0_dp) + (tzpx)*pp)*pp &
           &      -(tx*tz - 1.0_dp)*(0.5_dp*(tzmx**2 - 1.0_dp) + (tzmx)*mm)*mm &
           &    )/(2.0_dp*tx*tz)
        END IF
     ELSE
        IF(ABS(-bx-by+bz) < thr) THEN
           Wz = ( &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx+ty)*pp)*pp &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx-ty)*pm)*pm &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx+ty)*mp)*mp &
           &     +0.25_dp*(1.0_dp-ty)*(tx-1.0_dp)*(1.0_dp+tz)*(-ty) &
           &    )/(4.0_dp*tz*tx)
        ELSE IF(ABS(-bx+by+bz) < thr) THEN
           Wz = ( &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx+ty)*pp)*pp &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx-ty)*pm)*pm &
           &     +0.25_dp*(1.0_dp+ty)*(tx-1.0_dp)*(1.0_dp+tz)*ty &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx-ty)*mm)*mm &
           &    )/(4.0_dp*tz*tx)
        ELSE IF(ABS(bx-by+bz) < thr) THEN
           Wz = ( &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx+ty)*pp)*pp &
           &     +0.25_dp*(1.0_dp+ty)*(tx-1.0_dp)*(1.0_dp-tz)*ty &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx+ty)*mp)*mp &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx-ty)*mm)*mm &
           &    )/(4.0_dp*tz*tx)
        ELSE
           Wz = ( &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx+ty)*pp)*pp &
           &     -(tx*tz+1.0_dp)*(0.5_dp*(tzpx**2-1.0_dp)+(tzpx-ty)*pm)*pm &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx+ty)*mp)*mp &
           &     -(tx*tz-1.0_dp)*(0.5_dp*(tzmx**2-1.0_dp)+(tzmx-ty)*mm)*mm &
           &    )/(4.0_dp*tx*tz)
        END IF
     END IF
     !   
  ELSE IF(scdft_kernel == 2) THEN
     !
     txmy = TANH(0.5_dp*(x-y))
     !
     IF(bx < thr) THEN
        IF(ABS(by-bz) < thr) THEN
           Wz = (6.0_dp*tz-3.0_dp*z*(1.0_dp+tz**2)-z**3*(1.0_dp-tz**2))/(12.0_dp*z**3*tz)
        ELSE
           Wz = -( &
           &       ((( y+z)*(1.0_dp+ty*tz)-2.0_dp*( ty+tz)))*pp**3 &
           &      +(((-y+z)*(1.0_dp-ty*tz)-2.0_dp*(-ty+tz)))*mm**3 &
           &     ) /tz
        END IF
     ELSE IF(by < thr) THEN
        IF(ABS(bx-bz) < thr) THEN
           Wz = ( &
           &      -((tz + tx)*pp - 0.5_dp*(1.0_dp - tx**2))*pp &
           &      -0.25_dp*(1.0_dp - tx)*(1.0_dp + tx)*tx &
           &    )/(2.0_dp*tz*tx)
        ELSE
           Wz = ( &
           &      -((tz + tx)*pp - 0.5_dp*(1.0_dp - tx**2))*pp &
           &      +((tz - tx)*mm - 0.5_dp*(1.0_dp - tx**2))*mm &
           &    )/(2.0_dp*tz*tx)
        END IF
     ELSE
        IF(ABS(bx+by-bz) < thr) THEN
           Wz = ( &
           &      -(((1.0_dp+tx*ty)*tz+tx+ty)*pp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*pp &
           &      -(((1.0_dp-tx*ty)*tz+tx-ty)*pm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*pm &
           &      +(((1.0_dp-tx*ty)*tz-tx+ty)*mp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*mp &
           &      -0.25_dp*(1.0_dp-tx)*(1.0_dp-ty)*(1.0_dp+tx)*(1.0_dp+ty)*tx/(1.0_dp+tx*ty ) &
           &     )/(4.0_dp*tz*tx)
        ELSE IF(ABS(-bx+by+bz) < thr) THEN
           Wz = ( &
           &      -(((1.0_dp+tx*ty)*tz+tx+ty)*pp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*pp &
           &      -(((1.0_dp-tx*ty)*tz+tx-ty)*pm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*pm &
           &      -0.25_dp*(1.0_dp-tx)*(1.0_dp+ty)*tx*(1.0_dp+txmy) &
           &      +(((1.0_dp+tx*ty)*tz-tx-ty)*mm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*mm &
           &     )/(4.0_dp*tz*tx)
        ELSE IF(ABS(bx-by+bz) < thr) THEN
           Wz = ( &
           &      -(((1.0_dp+tx*ty)*tz+tx+ty)*pp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*pp &
           &      -0.25_dp*(1.0_dp-tx)*(1.0_dp+ty)*tx*(1.0_dp+txmy) &
           &      +(((1.0_dp-tx*ty)*tz-tx+ty)*mp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*mp &
           &      +(((1.0_dp+tx*ty)*tz-tx-ty)*mm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*mm &
           &     )/(4.0_dp*tz*tx)
        ELSE
           Wz = ( &
           &      -(((1.0_dp+tx*ty)*tz+tx+ty)*pp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*pp &
           &      -(((1.0_dp-tx*ty)*tz+tx-ty)*pm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*pm &
           &      +(((1.0_dp-tx*ty)*tz-tx+ty)*mp-0.5_dp*(1.0_dp-tx**2)*(1.0_dp+ty*tz))*mp &
           &      +(((1.0_dp+tx*ty)*tz-tx-ty)*mm-0.5_dp*(1.0_dp-tx**2)*(1.0_dp-ty*tz))*mm &
           &     )/(4.0_dp*tz*tx)
        END IF
     END IF
     ! 
  END IF
  !
  Wz = Wz * beta**2
  !
END FUNCTION Zweight
!
! Weight for K
!
FUNCTION Kweight(beta,x0,y0,z0) RESULT(Wk)
  !
  USE kinds, ONLY : DP
  USE sctk_val, ONLY : scdft_kernel, zero_kelvin
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: beta, x0, y0, z0
  REAL(dp) :: Wk
  !
  REAL(dp) :: thr = 1e-3_dp, tx, ty, tz, tg2z, x, y, z, bxymin, bx, by, bz
  REAL(dp) :: g1 = 1.33_dp, g2 = 3.88_dp
  !
  Wk = 0.0_dp
  !
  IF(zero_kelvin) THEN
     IF(scdft_kernel == 1) THEN
        Wk = - 2.0_dp / (x0+y0+z0)
     ELSE IF(scdft_kernel == 2) THEN
        Wk = -2d0*g1*g2*z0*(x0+y0+z0+g2*z0)/((x0+y0+z0)*(y0+g2*z0)*(x0+z0+g2*z0))
     END IF
     !
     RETURN
     !
  END IF
  !
  bx = beta * x0
  by = beta * y0
  bz = beta * z0
  !
  IF(scdft_kernel == 1) THEN
     !
     ! This kernel is symmetric for x and y
     !
     bxymin = MIN(bx, by)
     by = MAX(bx, by)
     bx = bxymin
     !
     IF(bx < thr) THEN
        IF(by < thr) THEN
           x = 0.0_dp
           y = 0.0_dp
           z = bz
        ELSE IF(ABS(by-bz) < thr) THEN
           x = 0.0_dp
           y = bz
           z = bz
        ELSE
           x = 0.0_dp
           y = by
           z = bz
        END IF
     ELSE
        IF(ABS(bx-by+bz) < thr) THEN
           x = bx
           y = bx + bz
           z = bz
        ELSE IF(ABS(bx+by-bz) < thr) THEN
           x = bx
           y = by
           z = bx + by
        ELSE
           x = bx
           y = by
           z = bz
        END IF
     END IF
     !
  ELSE IF(scdft_kernel == 2) THEN
     !
     IF(ABS(bx) < thr) THEN
        IF(ABS(by) < thr) THEN
           x = 0.0_dp
           y = 0.0_dp
           z = bz
        ELSE IF(ABS(by-bz) < thr) THEN
           x = 0.0_dp
           y = bz
           z = bz
        ELSE IF(ABS(by-g2*bz) < thr) THEN
           x = 0.0_dp
           y = g2*bz
           z = bz
        ELSE
           x = 0.0_dp
           y = by
           z = bz
        END IF
     ELSE IF(ABS(by) < thr) THEN
        IF(ABS(bx-bz) < thr) THEN
           x = bz
           y = 0.0_dp
           z = bz
        ELSE IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN
           x = (g2-1.0_dp)*bz
           y = 0.0_dp
           z = bz
        ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN
           x = (g2+1.0_dp)*bz
           y = 0.0_dp
           z = bz
        ELSE
           x = bx
           y = 0.0_dp
           z = bz
        END IF
     ELSE
        IF(ABS(by-g2*bz) < thr) THEN
           IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN ! x-y+z=0
              x = (g2-1.0_dp)*bz
              y = g2*bz
              z = bz
           ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN ! -x+y+z=0
              x = (g2+1.0_dp)*bz
              y = g2*bz
              z = bz
           ELSE
              x = bx
              y = g2*bz
              z = bz
           END IF 
        ELSE IF(ABS(bx+by-bz) < thr) THEN
           x = bx
           y = by
           z = bx + by
        ELSE IF(ABS(bx-by+bz) < thr) THEN
           x = bx
           y = bx + bz
           z = bz
        ELSE IF(ABS(-bx+by+bz) < thr) THEN
           x = by + bz
           y = by
           z = bz
        ELSE IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN
           x = (g2-1.0_dp)*bz
           y = by
           z = bz
        ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN
           x = (g2+1.0_dp)*bz
           y = by
           z = bz
        ELSE
           x = bx
           y = by
           z = bz
        END IF
     END IF
     !
  END IF
  !
  tx = TANH(0.5_dp*x)
  ty = TANH(0.5_dp*y)
  tz = TANH(0.5_dp*z)
  !
  IF(scdft_kernel == 1) THEN
     !
     IF(ABS(bx) < thr) THEN
        IF(ABS(by) < thr) THEN
           Wk = -2.0_dp*(4.0_dp*(2.0_dp/z-1.0_dp/tz)/z+1.0_dp)/z
        ELSE IF(ABS(by-bz) < thr) THEN
           Wk = -(1.0_dp+tz**2-2.0_dp*tz/z)/(2.0_dp*z*tz**2)
        ELSE
           Wk = ( &
           &     -(1.0_dp+ty*tz-2.0_dp*(tz+ty)/( y+z))/( y+z) &
           &     +(1.0_dp-ty*tz-2.0_dp*(tz-ty)/(-y+z))/(-y+z) &
           &    )/(ty*tz)
        END IF
     ELSE
        IF(ABS(bx-by+bz) < thr) THEN
           Wk = ( &
           &     -((1.0_dp+tx*ty)*tz+tx+ty)/( x+y+z) &
           &     +0.5_dp*(1.0_dp+tx)*(1.0_dp-ty)*(1.0_dp+tz) &
           &     +((1.0_dp-tx*ty)*tz-tx+ty)/(-x+y+z) &
           &     -((1.0_dp+tx*ty)*tz-tx-ty)/(-x-y+z) &
           &    )/(2.0_dp*tx*ty*tz)
        ELSE IF(ABS(bx+by-bz) < thr) THEN
           Wk = ( &
           &     -((1.0_dp+tx*ty)*tz+tx+ty)/( x+y+z) &
           &     +((1.0_dp-tx*ty)*tz+tx-ty)/( x-y+z) &
           &     +((1.0_dp-tx*ty)*tz-tx+ty)/(-x+y+z) &
           &     -0.5_dp*(1.0_dp-tx)*(1.0_dp-ty)*(1.0_dp+tz) &
           &    )/(2.0_dp*tx*ty*tz)
        ELSE
           Wk = ( &
           &     -((1.0_dp+tx*ty)*tz+tx+ty)/( x+y+z) &
           &     +((1.0_dp-tx*ty)*tz+tx-ty)/( x-y+z) &
           &     +((1.0_dp-tx*ty)*tz-tx+ty)/(-x+y+z) &
           &     -((1.0_dp+tx*ty)*tz-tx-ty)/(-x-y+z) &
           &    )/(2.0_dp*tx*ty*tz)
        END IF
     END IF
     !
  ELSE IF(scdft_kernel == 2) THEN
     !
     tg2z = TANH(0.5_dp*g2*z)
     !
     IF(ABS(bx) < thr) THEN
        IF(ABS(by) < thr) THEN
           Wk = 2.0_dp*g1* (-8.0_dp*g2**3*(g2**2-2.0_dp)-g2*(g2**2-1.0_dp)**2*z**2+2.0_dp*g2*(g2**2-1.0_dp)* &
           &    (2.0_dp *g2**2-1.0_dp)*z/tz &
           &    -2.0_dp*((g2**2-1.0_dp)* z+2.0_dp*(1.0_dp+g2**2)/tz)*tg2z)/(g2*(g2**2-1.0_dp)**2*z**3)
        ELSE IF(ABS(by-bz) < thr) THEN
           Wk = (g1*(-8.0_dp*(g2+g2**3)*tg2z+2.0_dp*g2**2*(9.0_dp-2.0_dp*g2**2+g2**4)*tz-g2*(g2**2-1.0_dp)* &
           &     (4.0_dp*tg2z*tz+g2**3*(1.0_dp+tz**2) &
           &     -g2*(5.0_dp+tz**2))*z))/(2.0_dp*(g2**2-1.0_dp)**3*tz**2*z**2)
        ELSE IF(ABS(by-g2*bz) < thr) THEN
           Wk = (g1*(4.0_dp*(3.0_dp*g2**4+6.0_dp*g2**2-1.0_dp)*tg2z-32.0_dp*g2**3*tz+2.0_dp*(g2**2-1.0_dp)* &
           &     (g2*(tg2z**2-1.0_dp+g2**2*(tg2z**2-3.0_dp)) &
           &    +(3.0_dp*g2**2-1.0_dp)*tg2z*tz)*z+g2*(g2**2-1.0_dp)**2*(tg2z**2-1.0_dp)*tz*z**2))/ &
           &    (2.0_dp*(g2**2-1.0_dp)**3*tg2z*tz*z**2)
        ELSE
           Wk = ( &
           &     -(-2.0_dp*( g2*z-y)*( y+(2.0_dp+g2)*z)+( y+z)*(( g2*z-y)*(1.0_dp+g2)*z/tz &
           &      -tg2z*( y+z)*((1.0_dp+g2)*z-2.0_dp/tz))+(1.0_dp+g2)**2*z**2*( y+z-2/tz)*ty) &
           &      /((g2*z-y)*(g2+1.0_dp)**2*( y+z)**2) &
           &     +(-2.0_dp*( g2*z+y)*(-y+(2.0_dp+g2)*z)+(-y+z)*(( g2*z+y)*(1.0_dp+g2)*z/tz &
           &      -tg2z*(-y+z)*((1.0_dp+g2)*z-2.0_dp/tz))+(1.0_dp+g2)**2*z**2*(-y+z-2/tz)*(-ty)) &
           &      /((g2*z+y)*(g2+1.0_dp)**2*(-y+z)**2) &
           &     -(-2.0_dp*(-g2*z-y)*( y+(2.0_dp-g2)*z)+( y+z)*((-g2*z-y)*(1.0_dp-g2)*z/tz &
           &      +tg2z*( y+z)*((1.0_dp-g2)*z-2.0_dp/tz))+(1.0_dp-g2)**2*z**2*( y+z-2/tz)*ty) &
           &      /((g2*z+y)*(g2-1.0_dp)**2*( y+z)**2) &
           &     +(-2.0_dp*(-g2*z+y)*(-y+(2.0_dp-g2)*z)+(-y+z)*((-g2*z+y)*(1.0_dp-g2)*z/tz &
           &      +tg2z*(-y+z)*((1.0_dp-g2)*z-2.0_dp/tz))+(1.0_dp-g2)**2*z**2*(-y+z-2/tz)*(-ty)) &
           &      /((g2*z-y)*(g2-1.0_dp)**2*(-y+z)**2) &
           &     )*g1*g2/(2.0_dp*ty*z)
        END IF
     ELSE IF(ABS(by) < thr) THEN
        IF(ABS(bx-bz) < thr) THEN
           Wk = (g1*(2.0_dp*g2**3*tz-8.0_dp*tg2z*(1.0_dp+tz**2)-g2*(g2**2-4.0_dp)*(1.0_dp+tz**2)*z)) &
           &     /(2.0_dp*g2*(g2**2-4.0_dp)*tz**2*z**2)
        ELSE IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN
           Wk = (g1*(-24.0_dp*tg2z*tz &
           &         + 2.0_dp*g2*(g2**3 + 20.0_dp*tg2z*tz - 9.0_dp*g2*tg2z*tz) &
           &         + g2*(g2**3* tg2z + (12.0_dp - 20.0_dp*g2 + 9.0_dp*g2**2 &
           &         + (-4.0_dp + (8.0_dp - 5.0_dp*g2)*g2)*tg2z**2)*tz - g2**3*tg2z*tz**2)*z &
           &         + (-2.0_dp*(-1.0_dp + g2)**2*(12.0_dp + (-4.0_dp + g2)*g2)*tg2z &
           &         + g2*(-4.0_dp + (8.0_dp - 5.0_dp*g2)*g2)*tg2z**2*z &
           &         + g2*(12.0_dp + g2*(-28.0_dp + g2*(21.0_dp + g2*(-5.0_dp + tz**2))))*z)/ tx)) &
           &  / (2.0_dp*(g2 - 2.0_dp)**2*(g2 - 1.0_dp)*g2**2*tz*z**2)
        ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN
           Wk = (2.0_dp*g1*(-g2**4*tz+g2*(2.0_dp+g2)**3*tg2z**2*tz+tg2z*(12.0_dp+28.0_dp*g2+21*g2**2+6.0_dp*g2**3 &
           &    -(12.0_dp+g2*(20+9*g2))*tz**2)) &
           &    -g1*g2*(1.0_dp+g2)*(2.0_dp+g2)*(2.0_dp*(tg2z**2-3.0_dp)*(tz**2-1.0_dp) &
           &    +g2*(5.0_dp+4.0_dp*tg2z*tz-tz**2+tg2z**2*(tz**2-1.0_dp)))*z) &
           &    /(2.0_dp*g2**2*(1.0_dp+g2)*(2.0_dp+g2)**2*tz*(tg2z+tz)*z**2)
        ELSE
           Wk = ( &
           &     (-2.0_dp*g2**2*z**2+(x-z)*( g2*z*( g2*z-x+z)+2.0_dp*(x-z)*tg2z   )*tx &
           &      +((x-z)*(  -g2*z *( g2*z-x+z)+2.0_dp*(-x+z)*tg2z   ) &
           &      +2.0_dp*g2**2*z**2*tx)/tz   )/(2.0_dp*g2**2*z**2*(x-z)**2*( g2*z-x+z)) &
           &    -(-2.0_dp*g2**2*z**2+(x+z)*( g2*z*( g2*z-x-z)+2.0_dp*(x+z)*tg2z   )*tx &
           &      +((x+z)*(  -g2*z *( g2*z-x-z)+2.0_dp*(-x-z)*tg2z   ) &
           &      +2.0_dp*g2**2*z**2*tx)/(-tz))/(2.0_dp*g2**2*z**2*(x+z)**2*( g2*z-x-z)) &
           &    -(-2.0_dp*g2**2*z**2+(x-z)*(-g2*z*(-g2*z-x+z)+2.0_dp*(x-z)*(-tg2z))*tx &
           &      +((x-z)*(-(-g2*z)*(-g2*z-x+z)+2.0_dp*(-x+z)*(-tg2z)) &
           &      +2.0_dp*g2**2*z**2*tx)/tz   )/(2.0_dp*g2**2*z**2*(x-z)**2*(-g2*z-x+z)) &
           &    +(-2.0_dp*g2**2*z**2+(x+z)*(-g2*z*(-g2*z-x-z)+2.0_dp*(x+z)*(-tg2z))*tx &
           &     +((x+z)*(-(-g2*z)*(-g2*z-x-z)+2.0_dp*(-x-z)*(-tg2z)) &
           &     +2.0_dp*g2**2*z**2*tx)/(-tz))/(2.0_dp*g2**2*z**2*(x+z)**2*(-g2*z-x-z)) &
           &    )*g1*g2*z/tx
        END IF
     ELSE
        IF(ABS(by-g2*bz) < thr) THEN
           IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN ! x-y+z=0
              Wk = ( &
              &     +(-2.0_dp*(tg2z+tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz+(tg2z**2-1.0_dp)*(tx*tz-1.0_dp)*(x+(g2-1.0_dp)*z)) &
              &      /(2.0_dp*(x+(g2-1.0_dp)*z)**2) &
              &     +((x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*(tx-ty-tz+tx*ty*tz)+y*(tx-tz+tg2z-tx*tz*tg2z)) &
              &      /((x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
              &     +tg2z*(1.0_dp-tg2z**2)*(1.0_dp-tz**2)/( 4.0_dp*(1-tz*tg2z)) &
              &     -(0.5_dp*(1.0_dp+tz*tx)*((-y-g2*z)*(1.0_dp-tg2z**2)+2.0_dp*tg2z)+ty*(1.0_dp+tz*tx)) &
              &      /((-x-y-z)*(-y-g2*z)) &
              &     -((x-z)*(tg2z+ty)*(1.0_dp-tx*tz)-g2*z*(tx+ty-tz-tx*ty*tz)-y*(tx-tz-tg2z+tx*tz*tg2z)) &
              &      /((x+y-z)*(x-(g2+1.0_dp)*z)*(y+g2*z)) &
              &     +(-2.0_dp*(-tg2z+tx)+2.0_dp*(1.0_dp-tg2z*tx)*tz+(tg2z**2-1.0_dp)*(tx*tz-1.0_dp)*(x+(-g2-1.0_dp)*z)) &
              &      /(2.0_dp*(x+(-g2-1.0_dp)*z)**2) &
              &     -(-2.0_dp*ty-2.0_dp*g2*(-tx+ty-tz)-tg2z*(2.0_dp+2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(-y-g2*(-x+y)) &
              &      -tx*ty*(-2.0_dp*(-1.0_dp-g2)*tz+(y+g2*(-x+y))+tz**2*(-y-g2*(-x+y)))) &
              &      /(2.0_dp*(-x-(g2+1.0_dp)*z)*(y+g2*z)) &
              &     -(-2.0_dp*(-tg2z-tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz+(tg2z**2-1.0_dp)*(-tx*tz-1.0_dp)*(-x+(-g2-1.0_dp)*z)) &
              &      /(2.0_dp*(-x+(-g2-1.0_dp)*z)**2) &
              &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
           ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN ! -x+y+z=0
              Wk = ( &
              &     +(-2.0_dp*(tg2z+tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz+(tg2z**2-1.0_dp)*(tx*tz-1.0_dp)*(x+(g2-1.0_dp)*z)) &
              &      /(2.0_dp*(x+(g2-1.0_dp)*z)**2) &
              &     -(2.0_dp*ty+2.0_dp*g2*(tx-ty-tz)+tg2z*(2.0_dp-2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(y+g2*(x-y)) &
              &      -tx*ty*(-2.0_dp*(-1.0_dp+g2)*tz+(-y-g2*(x-y))+tz**2*(y+g2*(x-y)))) &
              &      /(2.0_dp*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
              &     -(-2.0_dp*(tg2z-tx)+2.0_dp*(1.0_dp-tg2z*tx)*tz+(tg2z**2-1.0_dp)*(-tx*tz-1.0_dp)*(-x+(g2-1.0_dp)*z)) &
              &      /(2.0_dp*(-x+(g2-1.0_dp)*z)**2) &
              &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
              &      /((-x-y-z)*(x-(g2-1.0_dp)*z)*(-y-g2*z)) &
              &     +(0.5_dp*(1.0_dp-tx*tz)*((y+g2*z)*(1.0_dp-tg2z**2)-2.0_dp*tg2z)+ty*(tz*tx-1.0_dp))/((x+y-z)*(y+g2*z)) &
              &     +tg2z*(1.0_dp-tg2z**2)*(1.0_dp-tz**2)/( 4.0_dp*(1+tz*tg2z)) &
              &     +((-x-z)*(tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
              &     /((-x+y-z)*(-x-(g2+1.0_dp)*z)*(y+g2*z)) &
              &     -(-2.0_dp*(-tg2z-tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz &
              &     +(tg2z**2-1.0_dp)*(-tx*tz-1.0_dp)*(-x+(-g2-1.0_dp)*z))/(2.0_dp*(-x+(-g2-1.0_dp)*z)**2) &
              &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
           ELSE
              Wk = ( &
              &     +(-2.0_dp*(tg2z+tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz+(tg2z**2-1.0_dp)*(tx*tz-1.0_dp)*(x+(g2-1.0_dp)*z)) &
              &      /(2.0_dp*(x+(g2-1.0_dp)*z)**2) &
              &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*(tx-ty-tz+tx*ty*tz)+y*(tx-tz+tg2z-tx*tz*tg2z)) &
              &      /((x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
              &     -(-2.0_dp*(tg2z-tx)+2.0_dp*(1.0_dp-tg2z*tx)*tz+(tg2z**2-1.0_dp)*(-tx*tz-1.0_dp)*(-x+(g2-1.0_dp)*z)) &
              &      /(2.0_dp*(-x+(g2-1.0_dp)*z)**2) &
              &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
              &      /((-x-y-z)*(x-(g2-1.0_dp)*z)*(-y-g2*z)) &
              &     -(( x-z)*( tg2z+ty)*(1.0_dp-tx*tz)-g2*z*(tx+ty-tz-tx*ty*tz)-y*(tx-tz-tg2z+tx*tz*tg2z)) &
              &      /((x+y-z)*(x-(g2+1.0_dp)*z)*(y+g2*z)) &
              &     +(-2.0_dp*(-tg2z+tx)+2.0_dp*(1.0_dp-tg2z*tx)*tz+(tg2z**2-1.0_dp)*(tx*tz-1.0_dp)*(x+(-g2-1.0_dp)*z)) &
              &      /(2.0_dp*(x+(-g2-1.0_dp)*z)**2) &
              &     +((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
              &      /((-x+y-z)*(-x-(g2+1.0_dp)*z)*(y+g2*z)) &
              &     -(-2.0_dp*(-tg2z-tx)+2.0_dp*(1.0_dp+tg2z*tx)*tz+(tg2z**2-1.0_dp)*(-tx*tz-1.0_dp)*(-x+(-g2-1.0_dp)*z)) &
              &     /(2.0_dp*(-x+(-g2-1.0_dp)*z)**2) &
              &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
           END IF 
        ELSE IF(ABS(bx+by-bz) < thr) THEN
           Wk = ( &
           &     +(-2.0_dp*ty+2.0_dp*g2*(tx+ty-tz)+tg2z*(2.0_dp-2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(-y+g2*(x+y)) &
           &      +tx*ty*(-2.0_dp*(-1.0_dp+g2)*tz &
           &      +(y-g2*(x+y))+tz**2*(-y+g2*(x+y))))/(2.0_dp*(-x-(g2-1.0_dp)*z)*(y-g2*z)) &
           &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +((-x-z)*(-tg2z+ty)*(1.0_dp+tx*tz)+g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x+y-z)*( x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x-y-z)*( x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(-2.0_dp*ty-2.0_dp*g2*(tx+ty-tz)-tg2z*(2.0_dp-2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(-y-g2*(x+y)) &
           &      +tx*ty*(-2.0_dp*(-1.0_dp-g2)*tz &
           &      +(y+g2*(x+y))+tz**2*(-y-g2*(x+y))))/(2.0_dp*(x-(g2+1.0_dp)*z)*(y+g2*z)) &
           &     +((x-z)*(tg2z-ty)*(1.0_dp-tx*tz)-g2*z*(tx-ty-tz+tx*ty*tz)+y*(tx-tz-tg2z+tx*tz*tg2z)) &
           &      /((x-y-z)*(x-(g2+1.0_dp)*z)*(-y+g2*z)) &
           &     +((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x+y-z)*(-x-(g2+1.0_dp)*z)*( y+g2*z)) &
           &     -((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*(-x-(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        ELSE IF(ABS(bx-by+bz) < thr) THEN
           Wk = ( &
           &     -(( x-z)*(-tg2z+ty)*(1.0_dp-tx*tz)+g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     -(-2.0_dp*ty+2.0_dp*g2*(-tx+ty-tz)+tg2z*(2.0_dp+2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(-y+g2*(-x+y)) &
           &      -tx*ty*(-2.0_dp*(-1.0_dp+g2)*tz+(y-g2*(-x+y))+tz**2*(-y+g2*(-x+y))))/(2.0_dp*(x-(g2-1.0_dp)*z)*(y-g2*z)) &
           &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x-y-z)*( x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(( x-z)*( tg2z+ty)*(1.0_dp-tx*tz)-g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     -(( x-z)*( tg2z-ty)*(1.0_dp-tx*tz)-g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &     -(-2.0_dp*ty-2.0_dp*g2*(-tx+ty-tz)-tg2z*(2.0_dp+2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(-y-g2*(-x+y)) &
           &      -tx*ty*(-2.0_dp*(-1.0_dp-g2)*tz+(y+g2*(-x+y))+tz**2*(-y-g2*(-x+y)))) &
           &      /(2.0_dp*(-x-(g2+1.0_dp)*z)*(y+g2*z)) &
           &     +((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*( x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        ELSE IF(ABS(-bx+by+bz) < thr) THEN
           Wk = ( &
           &     -(( x-z)*(-tg2z+ty)*(1.0_dp-tx*tz)+g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     -(2.0_dp*ty+2.0_dp*g2*(tx-ty-tz)+tg2z*(2.0_dp-2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(y+g2*(x-y)) &
           &      -tx*ty*(-2.0_dp*(-1.0_dp+g2)*tz+(-y-g2*(x-y))+tz**2*(y+g2*(x-y))))/(2.0_dp*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +((-x-z)*(-tg2z+ty)*(1.0_dp+tx*tz)+g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x+y-z)*( x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x-y-z)*( x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(( x-z)*( tg2z+ty)*(1.0_dp-tx*tz)-g2*z*(tx+ty-tz-tx*ty*tz)-y*(tx-tz-tg2z+tx*tz*tg2z)) &
           &      /((x+y-z)*(-x+(g2+1.0_dp)*z)*(y+g2*z)) &
           &     -(2.0_dp*ty-2.0_dp*g2*(tx-ty-tz)-tg2z*(2.0_dp-2.0_dp*tx*tz)+(-1.0_dp+tz**2)*(y-g2*(x-y)) &
           &      -tx*ty*(-2.0_dp*(-1.0_dp-g2)*tz+(-y+g2*(x-y))+tz**2*(y-g2*(x-y))))/(2.0_dp*(x-(g2+1.0_dp)*z)*(-y+g2*z)) &
           &     -((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x+y-z)*( x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     +((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*( x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        ELSE IF(ABS(bx-(g2-1.0_dp)*bz) < thr) THEN
           Wk = ( &
           &     -(( x-z)*(-tg2z+ty)*(1.0_dp-tx*tz)+g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(0.5_dp*(1.0_dp+tz*tx)*(( y-g2*z)*(1.0_dp-tg2z**2)+2.0_dp* tg2z)-ty*(1.0_dp+tz*tx))/((-x+y-z)*( y-g2*z)) &
           &     -(0.5_dp*(1.0_dp+tz*tx)*((-y-g2*z)*(1.0_dp-tg2z**2)+2.0_dp* tg2z)+ty*(1.0_dp+tz*tx))/((-x-y-z)*(-y-g2*z)) &
           &     +(( x-z)*( tg2z+ty)*(1.0_dp-tx*tz)-g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     -(( x-z)*( tg2z-ty)*(1.0_dp-tx*tz)-g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &     -((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x+y-z)*( x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     +((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*( x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        ELSE IF(ABS(bx-(g2+1.0_dp)*bz) < thr) THEN
           Wk = ( &
           &     -(( x-z)*(-tg2z+ty)*(1.0_dp-tx*tz)+g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &       /(( x+y-z)*(-x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +((-x-z)*(-tg2z+ty)*(1.0_dp+tx*tz)+g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x+y-z)*( x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x-y-z)*( x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(0.5_dp*(1.0_dp-tx*tz)*(( y+g2*z)*(1.0_dp-tg2z**2)-2.0_dp*tg2z)+ty*(tz*tx-1.0_dp)) &
           &      /((x+y-z)*( y+g2*z)) &
           &     -(0.5_dp*(1.0_dp-tx*tz)*((-y+g2*z)*(1.0_dp-tg2z**2)-2.0_dp*tg2z)-ty*(tz*tx-1.0_dp)) &
           &      /((x-y-z)*(-y+g2*z)) &
           &     -((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x+y-z)*( x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     +((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*( x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        ELSE
           Wk = ( &
           &     -(( x-z)*(-tg2z+ty)*(1.0_dp-tx*tz)+g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     +(( x-z)*(-tg2z-ty)*(1.0_dp-tx*tz)+g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz+tg2z-tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +((-x-z)*(-tg2z+ty)*(1.0_dp+tx*tz)+g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x+y-z)*( x-(g2-1.0_dp)*z)*( y-g2*z)) &
           &     -((-x-z)*(-tg2z-ty)*(1.0_dp+tx*tz)+g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz+tg2z+tx*tz*tg2z)) &
           &      /((-x-y-z)*( x-(g2-1.0_dp)*z)*(-y-g2*z)) &
           &     +(( x-z)*( tg2z+ty)*(1.0_dp-tx*tz)-g2*z*( tx+ty-tz-tx*ty*tz)-y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x+y-z)*(-x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     -(( x-z)*( tg2z-ty)*(1.0_dp-tx*tz)-g2*z*( tx-ty-tz+tx*ty*tz)+y*( tx-tz-tg2z+tx*tz*tg2z)) &
           &      /(( x-y-z)*(-x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &     -((-x-z)*( tg2z+ty)*(1.0_dp+tx*tz)-g2*z*(-tx+ty-tz+tx*ty*tz)-y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x+y-z)*( x+(g2+1.0_dp)*z)*( y+g2*z)) &
           &     +((-x-z)*( tg2z-ty)*(1.0_dp+tx*tz)-g2*z*(-tx-ty-tz-tx*ty*tz)+y*(-tx-tz-tg2z-tx*tz*tg2z)) &
           &      /((-x-y-z)*( x+(g2+1.0_dp)*z)*(-y+g2*z)) &
           &    )*g1*g2*z/(4.0_dp*tx*ty*tz)
        END IF
     END IF
     !
  END IF
  !
  Wk = Wk * beta
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
END MODULE sctk_kernel_weight
