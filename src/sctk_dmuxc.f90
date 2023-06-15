!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_dmuxc
  !
  IMPLICIT NONE
  !
CONTAINS
!
!--------------------------------------------------------
SUBROUTINE generate_dmuxc()
  !------------------------------------------------------
  !
  ! This routine output the f_{XC} for LDA in G space to a file.
  !
  USE kinds, ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE eqv,       ONLY : dmuxc
  USE scf,              ONLY : rho, rho_core
  USE noncollin_module, ONLY : nspin_gga, nspin_mag, noncolin
  USE gvect, ONLY : ngm
  USE spin_orb, ONLY : domag
  !
  REAL(DP) :: rho_of_r(dfftp%nnr)
  COMPLEX(DP) :: rho_of_g(ngm)
  !
  ! Computes the derivative of the XC potential
  !
  nspin_gga = 2
  domag = .TRUE.
  nspin_mag = 4
  noncolin = .TRUE.
  IF(.NOT. ALLOCATED(dmuxc)) ALLOCATE(dmuxc(dfftp%nnr, nspin_mag, nspin_mag))
  !
  rho_of_r(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr, 1)
  rho_of_g(1:ngm) = rho%of_g(1:ngm, 1)
  !
  DEALLOCATE(rho%of_r, rho%of_g)
  ALLOCATE(rho%of_r(dfftp%nnr, nspin_mag), rho%of_g(ngm,nspin_mag))
  !
  rho%of_r(1:dfftp%nnr, 1) = rho_of_r(1:dfftp%nnr)
  rho%of_r(1:dfftp%nnr, 2:4) = 0.0_dp
  !
  rho%of_g(1:ngm, 1) = rho_of_g(1:ngm)
  rho%of_g(1:ngm, 2:nspin_mag) = 0.0_dp
  !
  ! LDA
  !
  rho_of_r(1:dfftp%nnr) = rho_of_r(1:dfftp%nnr) + rho_core(1:dfftp%nnr)
  !
  CALL dmxc_nc_para(dfftp%nnr, rho_of_r, dmuxc)
  !
  ! GGA
  !
  CALL setup_dgc()
  !CALL compute_vsgga2()
  !
END SUBROUTINE generate_dmuxc
!
! Apply XC term
!
SUBROUTINE apply_xc()
  !
  USE kinds, ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE eqv,               ONLY : dmuxc
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE funct,             ONLY : dft_is_gradient
  USE noncollin_module,  ONLY : nspin_mag
  !
  USE sctk_invert, ONLY : invert
  USE sctk_val, ONLY : gindx, gq2, laddxc, nf, ngv, ngv0, ngv1, nmf, wscr
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, imf, igv(3), gindx_p(ngv)
  COMPLEX(dp) :: vec(dfftp%nnr), drho(dfftp%nnr, nspin_mag), dvgga(dfftp%nnr, nspin_mag)
  !
  IF(laddxc>0) THEN
     !
     CALL start_clock("apply_xc")
     !
     DO ig = 1, ngv
        !
        igv(3) = (gindx(ig) - 1) / (nf(1)*nf(2))
        igv(2) = (gindx(ig) - 1 - igv(3)*nf(2)*nf(1)) / nf(1)
        igv(1) =  gindx(ig) - 1 - igv(3)*nf(2)*nf(1) - igv(2)*nf(1)
        WHERE(igv(1:3)*2 >= nf(1:3)) igv(1:3) = igv(1:3) - nf(1:3)
        !
        igv(1:3) = modulo(igv(1:3), (/dfftp%nr1, dfftp%nr2, dfftp%nr3/))
        gindx_p(ig) = 1 + igv(1) + igv(2) * dfftp%nr1x + igv(3) * dfftp%nr1x * dfftp%nr2x
        !
     END DO
     !
     DO imf = 0, nmf
        DO ig = ngv0, ngv1
           !
           vec(1:dfftp%nnr) = 0.0_dp
           vec(gindx_p(1:ngv)) = wscr(1:ngv, ig, imf, 1)
           !
           ! LDA/GGA kernel is computed in the real space
           !
           CALL invfft ('Rho', vec, dfftp)
           drho(1:dfftp%nnr, 1) = vec(1:dfftp%nnr)
           drho(1:dfftp%nnr, 2:nspin_mag) = 0.0_dp
           !
           vec(1:dfftp%nnr) = vec(1:dfftp%nnr) * dmuxc(1:dfftp%nnr,1,1)
           !
           IF(dft_is_gradient()) THEN
              dvgga(1:dfftp%nnr, 1:4) = 0.0_dp
              CALL dgradcorr2(drho, dvgga)
              vec(1:dfftp%nnr) = vec(1:dfftp%nnr) + dvgga(1:dfftp%nnr, 1)
           END IF
           !              
           CALL fwfft ('Rho', vec, dfftp)
           !
           wscr(1:ngv, ig, imf, 1) = wscr(1:ngv, ig, imf, 1) + gq2(1:ngv) * vec(gindx_p(1:ngv))
           !
        END DO ! ig
     END DO ! imf
     !
     CALL stop_clock("apply_xc")
     !
  END IF
  !
  wscr(1:ngv, ngv0:ngv1, 0:nmf, 1) = - wscr(1:ngv, ngv0:ngv1, 0:nmf, 1)
  DO ig = ngv0, ngv1
     wscr(ig, ig, 0:nmf, 1) = CMPLX(gq2(ig), 0.0_dp, KIND=dp) + wscr(ig, ig, 0:nmf, 1)
  END DO
  CALL invert(wscr(1:ngv, ngv0:ngv1, 0:nmf, 1))
  !
END SUBROUTINE apply_xc
!
!
!
SUBROUTINE apply_xc_spin()
  !
  USE kinds, ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE eqv,               ONLY : dmuxc
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE funct,             ONLY : dft_is_gradient
  USE noncollin_module,  ONLY : nspin_mag, npol
  !
  USE sctk_val, ONLY : gindx, nf, nmf, ngv, ngv0, ngv1, wscr
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, imf, igv(3), gindx_p(ngv), ipol
  COMPLEX(dp) :: vec(dfftp%nnr), drho(dfftp%nnr, nspin_mag)!, dvgga(dfftp%nnr, nspin_mag)
  !
  DO ig = 1, ngv
     !
     igv(3) = (gindx(ig) - 1) / (nf(1)*nf(2))
     igv(2) = (gindx(ig) - 1 - igv(3)*nf(2)*nf(1)) / nf(1)
     igv(1) =  gindx(ig) - 1 - igv(3)*nf(2)*nf(1) - igv(2)*nf(1)
     WHERE(igv(1:3)*2 >= nf(1:3)) igv(1:3) = igv(1:3) - nf(1:3)
     !
     igv(1:3) = modulo(igv(1:3), (/dfftp%nr1, dfftp%nr2, dfftp%nr3/))
     gindx_p(ig) = 1 + igv(1) + igv(2) * dfftp%nr1x + igv(3) * dfftp%nr1x * dfftp%nr2x
     !
  END DO
  !
  DO ipol = 2, 2*npol
     !
     DO imf = 0, nmf
        !
        DO ig = ngv0, ngv1
           !
           vec(1:dfftp%nnr) = 0.0_dp
           vec(gindx_p(1:ngv)) = wscr(1:ngv, ig, imf, ipol)
           !
           ! LDA/GGA kernel is computed in the real space
           !
           CALL invfft ('Rho', vec, dfftp)
           drho(1:dfftp%nnr, 1:4) = 0.0_dp
           drho(1:dfftp%nnr, ipol) = vec(1:dfftp%nnr)
           !
           vec(1:dfftp%nnr) = drho( 1:dfftp%nnr, ipol) &
           &                * dmuxc(1:dfftp%nnr, ipol, ipol)
           !
           !IF(dft_is_gradient()) THEN ! tentative: sGGA turn off
           !   CALL dgradcorr2(drho, dvgga) 
           !   vec(1:dfftp%nnr) = vec(1:dfftp%nnr) + dvgga(1:dfftp%nnr, ipol)
           !END IF
           !              
           CALL fwfft ('Rho', vec, dfftp)
           !
           wscr(1:ngv, ig, imf, ipol) = vec(gindx_p(1:ngv))
           !
        END DO ! ig
        !
     END DO ! ipol
     !
  END DO ! imf = 0, nmf
  !
END SUBROUTINE apply_xc_spin
!
!
!
SUBROUTINE dmxc_nc_para( length, rho_in, dmuxc )
!-----------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density and magnetization in the non-collinear case.
  !
  USE xc_lda_lsda,  ONLY: xc_lsda
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
  !! magnetization vector
  REAL(DP), INTENT(OUT), DIMENSION(length,4,4) :: dmuxc
  !! derivative of XC functional
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length)     :: rhotot, dr, dz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rhoaux, zetaux
  REAL(DP), DIMENSION(length) :: vs
  INTEGER, DIMENSION(length) :: null_v
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: aux1, aux2
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vx, vc
  REAL(DP), DIMENSION(length) :: dvxc_rho, dbx_rho, dby_rho, dbz_rho
  !
  REAL(DP) :: dvxc_mx, dvxc_my, dvxc_mz, &
  &           dbx_mx, dbx_my, dbx_mz,    &
  &           dby_mx, dby_my, dby_mz,    &
  &           dbz_mx, dbz_my, dbz_mz
  !
  INTEGER :: i1, i2, i3, i4, i5, i
  INTEGER :: f1, f2, f3, f4, f5
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP, &
  &                      rho_trash = 0.5_DP
  !
  dmuxc = 0.0_DP
  !
  ALLOCATE( rhoaux(length*5), zetaux(length*5) )
  ALLOCATE( aux1(length*5), aux2(length*5) )
  ALLOCATE( vx(length*5,2), vc(length*5,2) )
  !
  rhotot = rho_in
  null_v = 1
  !
  i1 = 1     ;   f1 = length    !           five blocks:  [ rho    , zeta    ]   
  i2 = f1+1  ;   f2 = 2*length  !                         [ rho+dr , zeta    ]   
  i3 = f2+1  ;   f3 = 3*length  !                         [ rho-dr , zeta    ]   
  i4 = f3+1  ;   f4 = 4*length  !                         [ rho    , zeta+dz ]   
  i5 = f4+1  ;   f5 = 5*length  !                         [ rho    , zeta-dz ]   
  !
  dz = 1.0E-6_DP     !dz = MIN( 1.d-6, 1.d-4*ABS(zeta) ) 
  !
  DO i = 1, length
     IF (rhotot(i) <= small) THEN
        rhotot(i) = rho_trash ; null_v(i) = 0.0_DP
     ENDIF
  ENDDO
  !
  dr = MIN( 1.E-6_DP, 1.E-4_DP * rhotot )
  !   
  rhoaux(i1:f1) = rhotot         ;   zetaux(i1:f1) = 0.0_DP
  rhoaux(i2:f2) = rhotot + dr    ;   zetaux(i2:f2) = 0.0_DP
  rhoaux(i3:f3) = rhotot - dr    ;   zetaux(i3:f3) = 0.0_DP
  rhoaux(i4:f4) = rhotot         ;   zetaux(i4:f4) = + dz
  rhoaux(i5:f5) = rhotot         ;   zetaux(i5:f5) = - dz
  !
  !
  CALL xc_lsda( length*5, rhoaux, zetaux, aux1, aux2, vx, vc )
  !
  !
  vs(:) = 0.5_DP*( vx(i1:f1,1)+vc(i1:f1,1)-vx(i1:f1,2)-vc(i1:f1,2) )
  !
  dvxc_rho(:) = ((vx(i2:f2,1) + vc(i2:f2,1) - vx(i3:f3,1) - vc(i3:f3,1)) + &
  &              (vx(i2:f2,2) + vc(i2:f2,2) - vx(i3:f3,2) - vc(i3:f3,2))) / (4.0_DP*dr)
  !   
  aux2(1:length) =  vx(i2:f2,1) + vc(i2:f2,1) - vx(i3:f3,1) - vc(i3:f3,1) - &
  &               ( vx(i2:f2,2) + vc(i2:f2,2) - vx(i3:f3,2) - vc(i3:f3,2) )
  !
  dbx_rho(:) = aux2(1:length) / (4.0_DP*dr)
  dby_rho(:) = aux2(1:length) / (4.0_DP*dr)
  dbz_rho(:) = aux2(1:length) / (4.0_DP*dr)
  !
  aux1(1:length) =  vx(i4:f4,1) + vc(i4:f4,1) - vx(i5:f5,1) - vc(i5:f5,1) + &
  &                 vx(i4:f4,2) + vc(i4:f4,2) - vx(i5:f5,2) - vc(i5:f5,2)
  aux2(1:length) =  vx(i4:f4,1) + vc(i4:f4,1) - vx(i5:f5,1) - vc(i5:f5,1) - &
  &               ( vx(i4:f4,2) + vc(i4:f4,2) - vx(i5:f5,2) - vc(i5:f5,2) )
  !
  DO i = 1, length
     !
     IF ( null_v(i) == 0 ) THEN
        dmuxc(i,:,:) = 0.0_DP
        CYCLE
     ENDIF
     !
     dmuxc(i,1,1) = dvxc_rho(i)
     dmuxc(i,2,1) = dbx_rho(i)
     dmuxc(i,3,1) = dby_rho(i)
     dmuxc(i,4,1) = dbz_rho(i)
     !
     ! ... Here the derivatives with respect to m
     !
     dvxc_mx = aux1(i) / rhotot(i) / (4.0_DP*dz(i))
     dvxc_my = aux1(i) / rhotot(i) / (4.0_DP*dz(i))
     dvxc_mz = aux1(i) / rhotot(i) / (4.0_DP*dz(i))
     !
     dbx_mx  = aux2(i) / rhotot(i) / (4.0_DP*dz(i))
     dbx_my  = 0.0_DP
     dbx_mz  = 0.0_DP
     !
     dby_mx  = 0.0_DP
     dby_my  = aux2(i) / rhotot(i) / (4.0_DP*dz(i))
     dby_mz  = 0.0_DP
     !
     dbz_mx  = 0.0_DP
     dbz_my  = 0.0_DP
     dbz_mz  = aux2(i) / rhotot(i) / (4.0_DP*dz(i))
     !
     ! ... assigns values to dmuxc and sets to zero trash points
     !
     dmuxc(i,1,2) = dvxc_mx 
     dmuxc(i,1,3) = dvxc_my  
     dmuxc(i,1,4) = dvxc_mz 
     !
     dmuxc(i,2,2) = dbx_mx 
     dmuxc(i,2,3) = dbx_my 
     dmuxc(i,2,4) = dbx_mz 
     !
     dmuxc(i,3,2) = dby_mx 
     dmuxc(i,3,3) = dby_my 
     dmuxc(i,3,4) = dby_mz 
     !
     dmuxc(i,4,2) = dbz_mx 
     dmuxc(i,4,3) = dbz_my 
     dmuxc(i,4,4) = dbz_mz 
     !
  ENDDO
  !
  ! ... brings to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  DEALLOCATE( rhoaux, zetaux)
  DEALLOCATE( aux1, aux2 )
  DEALLOCATE( vx, vc )
  !
  RETURN
  !
END SUBROUTINE dmxc_nc_para
!
!
!
!!SUBROUTINE compute_vsgga2()
!!  !----------------------------------------------------------------------------
!!  !
!!  USE constants,            ONLY : e2
!!  USE kinds,                ONLY : DP
!!  USE gvect,                ONLY : g
!!  USE noncollin_module,     ONLY : nspin_gga
!!  USE funct,                ONLY : gcx_spin, gcc_spin, &
!!                                   dft_is_gradient, get_igcc
!!  USE fft_base,             ONLY : dfftp
!!  USE gc_lr,            ONLY : grho, vsgga
!!  USE scf, ONLY : rho
!!  !
!!  IMPLICIT NONE
!!  !
!!  INTEGER :: k, ipol, is
!!  !
!!  REAL(DP),    ALLOCATABLE :: h(:,:,:), dh(:)
!!  REAL(DP),    ALLOCATABLE :: vaux(:,:)
!!  !
!!  LOGICAL  :: igcc_is_lyp
!!  REAL(DP) :: grho2(2), sx, sc, v2c, &
!!              v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw , &
!!              arho, rh, grh2
!!  REAL(DP) :: v2cup, v2cdw, v2cud, grup, grdw, dz = 1.0E-6_DP
!!  !
!!  REAL(DP), PARAMETER :: vanishing_charge = 1.D-6
!!  REAL(DP), PARAMETER :: epsg = 1.D-10
!!  !
!!  igcc_is_lyp = (get_igcc() == 3)
!!  !
!!  ALLOCATE(    h( 3, dfftp%nnr, nspin_gga) )
!!  ALLOCATE( vaux( dfftp%nnr, nspin_gga ) )
!!  !
!!  DO k = 1, dfftp%nnr
!!     !
!!     rh = rho%of_r(k,1)
!!     arho = ABS(rh)
!!     !
!!     IF ( arho > vanishing_charge ) THEN
!!        !
!!        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
!!        !
!!        IF ( grho2(1) > epsg .OR. grho2(2) > epsg ) THEN
!!           CALL gcx_spin(rh*(1.0_dp+dz)*0.5_dp, rh*(1.0_dp-dz)*0.5_dp, &
!!           &             grho2(1), grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
!!           !
!!           v1xup = v1xup / (rh*dz)
!!           v1xdw = v1xdw / (rh*dz)
!!           v2xup = v2xup / (rh*dz)
!!           v2xdw = v2xdw / (rh*dz)
!!           !
!!           grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 &
!!           &    + ( grho(2,k,1) + grho(2,k,2) )**2 &
!!           &    + ( grho(3,k,1) + grho(3,k,2) )**2
!!           !
!!           CALL gcc_spin( rh, dz, grh2, sc, v1cup, v1cdw, v2c )
!!           !
!!           v1cup = v1cup / (rh*dz)
!!           v1cdw = v1cdw / (rh*dz)
!!           v2c = v2c / (rh*dz)
!!           v2cup = v2c
!!           v2cdw = v2c
!!           v2cud = v2c
!!           !
!!        ELSE
!!           !
!!           sc    = 0.D0
!!           sx    = 0.D0
!!           v1xup = 0.D0
!!           v1xdw = 0.D0
!!           v2xup = 0.D0
!!           v2xdw = 0.D0
!!           v1cup = 0.D0
!!           v1cdw = 0.D0
!!           v2c   = 0.D0
!!           v2cup = 0.D0
!!           v2cdw = 0.D0
!!           v2cud = 0.D0
!!        ENDIF
!!     ELSE
!!        !
!!        sc    = 0.D0
!!        sx    = 0.D0
!!        v1xup = 0.D0
!!        v1xdw = 0.D0
!!        v2xup = 0.D0
!!        v2xdw = 0.D0
!!        v1cup = 0.D0
!!        v1cdw = 0.D0
!!        v2c   = 0.D0
!!        v2cup = 0.D0
!!        v2cdw = 0.D0
!!        v2cud = 0.D0
!!        !
!!     ENDIF
!!     !
!!     ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
!!     !
!!     vaux(k,1) = e2 * ( v1xup + v1cup )
!!     vaux(k,2) = e2 * ( v1xdw + v1cdw )
!!     !
!!     ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
!!     !
!!     DO ipol = 1, 3
!!        !
!!        grup = grho(ipol,k,1)
!!        grdw = grho(ipol,k,2)
!!        h(ipol,k,1) = e2 * ( ( v2xup + v2cup ) * grup + v2cud * grdw )
!!        h(ipol,k,2) = e2 * ( ( v2xdw + v2cdw ) * grdw + v2cud * grup )
!!        !
!!     END DO
!!     !
!!  END DO
!!  !
!!  ALLOCATE( dh( dfftp%nnr ) )
!!  !
!!  ! ... second term of the gradient correction :
!!  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
!!  !
!!  DO is = 1, nspin_gga
!!     !
!!     CALL fft_graddot(dfftp, h(1,1,is), g, dh )
!!     !
!!     vaux(:,is) = vaux(:,is) - dh(:)
!!     !
!!  END DO
!!  !
!!  vsgga(:)=(vaux(:,1)-vaux(:,2))
!!  !
!!  DEALLOCATE( dh )
!!  DEALLOCATE( h )
!!  DEALLOCATE( vaux )
!!  !
!!  RETURN
!!  !
!!END SUBROUTINE compute_vsgga2
!
!
!
subroutine dgradcorr2 (drho, dvxc)
  !--------------------------------------------------------------------
  !
  !  Add gradient correction contribution to 
  !  the responce exchange-correlation potential dvxc.
  !
  USE fft_base,  ONLY : dfftp
  USE kinds,            ONLY : DP
  USE gc_lr,            ONLY : vsgga, segni
  USE fft_types,        ONLY : fft_type_descriptor
  USE gc_lr,            ONLY : segni, vsgga, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE noncollin_module, ONLY : nspin_gga, nspin_mag
  USE disp,              ONLY : x_q
  USE scf,               ONLY : rho
  USE control_ph, ONLY : current_iq
  USE gvect, ONLY : g
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(IN) :: drho(dfftp%nnr, nspin_mag)
  COMPLEX(DP),INTENT(OUT) :: dvxc(dfftp%nnr, nspin_mag)
  !
  REAL(DP) :: grho2, seg, seg0
  COMPLEX(DP) :: term
  COMPLEX(DP) :: a(2, 2, 2), b(2, 2, 2, 2), c(2, 2, 2), &
  &              ps(2, 2), ps1(3, 2, 2), ps2(3, 2, 2, 2)
  COMPLEX(DP),ALLOCATABLE :: gdrho(:,:,:), h(:,:,:), dh(:)
  COMPLEX(DP),ALLOCATABLE :: gdmag(:,:,:), vgg(:,:)
  COMPLEX(DP),ALLOCATABLE :: drhoout(:,:)
  REAL(DP),ALLOCATABLE :: rhoout(:,:)
  INTEGER :: k, ipol, is, js, ks, ls
  !
  ALLOCATE(gdmag(3, dfftp%nnr, nspin_mag))
  ALLOCATE(vgg(dfftp%nnr, nspin_gga))
  vgg = (0.0_dp,0.0_dp)
  ALLOCATE(rhoout(dfftp%nnr, nspin_gga))
  ALLOCATE(drhoout(dfftp%nnr, nspin_gga))
  ALLOCATE(gdrho(3, dfftp%nnr, nspin_gga))
  ALLOCATE(h(3, dfftp%nnr, nspin_gga))
  ALLOCATE(dh(dfftp%nnr))
  !
  h(:, :, :) = (0.d0, 0.d0)
  DO is = 1, nspin_mag
     CALL fft_qgradient(dfftp, drho(1,is), x_q(1:3,current_iq), g, gdmag(1, 1, is) )
  END DO
  DO is = 1, nspin_gga
     IF(is==1) seg0 =  0.5_dp
     IF(is==2) seg0 = -0.5_dp
     rhoout( :,is) = 0.5_dp*rho%of_r( :,1)
     drhoout(:,is) = 0.5_dp*drho(:,1)
     DO ipol = 1, 3
        gdrho(ipol,:,is) = 0.5_dp*gdmag(ipol,:,1)
     END DO
     DO k = 1, dfftp%nnr
        seg = seg0*segni(k)
        drhoout(k,is) = drhoout(k,is) &
        &             + seg*SQRT(DOT_PRODUCT(drho(k,2:4), drho(k,2:4)))
        DO ipol = 1,3
           gdrho(ipol,k,is) = gdrho(ipol,k,is) &
           &                + seg*SQRT(DOT_PRODUCT(gdmag(ipol,k,2:4), gdmag(ipol,k,2:4)))
        END DO
     END DO
  END DO
  !
  DO k = 1, dfftp%nnr
     grho2 = grho(1, k, 1)**2 + grho(2, k, 1)**2 + grho(3, k, 1)**2
     !
     ps(:,:) = (0.d0, 0.d0)
     DO is = 1, nspin_gga
        DO js = 1, nspin_gga
           ps1(1:3, is, js) = drhoout(k, is) * grho(1:3, k, js)
           ps(is, js) = ps(is, js) + SUM(grho(1:3,k,is)*gdrho(1:3,k,js))
           DO ks = 1, nspin_gga
              IF(is == js .AND. js == ks) THEN
                 a(is, js, ks) = dvxc_sr(k, is, is)
                 c(is, js, ks) = dvxc_sr(k, is, is)
              ELSE
                 IF(is == 1) THEN
                    a(is, js, ks) = dvxc_sr(k, 1, 2)
                 ELSE
                    a(is, js, ks) = dvxc_sr(k, 2, 1)
                 END IF
                 IF(js == 1) THEN
                    c(is, js, ks) = dvxc_sr(k, 1, 2)
                 ELSE
                    c(is, js, ks) = dvxc_sr(k, 2, 1)
                 END IF
              END IF
              ps2(1:3, is, js, ks) = ps(is, js) * grho(1:3, k, ks)
              DO ls = 1, nspin_gga
                 IF(is == js .AND. js == ks .AND. ks == ls) THEN
                    b(is, js, ks, ls) = dvxc_ss(k, is, is)
                 ELSE
                    IF(is == 1) THEN
                       b(is, js, ks, ls) = dvxc_ss(k, 1, 2)
                    ELSE
                       b(is, js, ks, ls) = dvxc_ss(k, 2, 1)
                    END IF
                 END IF
              END DO
           END DO
        END DO
     END DO
     DO is = 1, nspin_gga
        DO js = 1, nspin_gga
           vgg(k, is) = vgg(k, is) + dvxc_rr(k,is,js)*drhoout(k, js)
           h(1:3, k, is) = h(1:3, k, is) &
           &             + dvxc_s(k, is, js) * gdrho(1:3, k, js)
           DO ks = 1, nspin_gga
              vgg(k, is) = vgg(k, is) + a(is, js, ks) * ps(js, ks)
              h(1:3, k, is) = h(1:3, k, is) &
              &             + c(is, js, ks) * ps1(1:3, js, ks)
              DO ls = 1, nspin_gga
                 h(1:3, k, is) = h(1:3, k, is) &
                 &             + b(is, js, ks, ls) * ps2(1:3, js, ks, ls)
              END DO
           END DO
        END DO
     END DO
  END DO
  ! linear variation of the second term
  DO is = 1, nspin_gga
     CALL fft_qgraddot(dfftp, h(1, 1, is), x_q(1:3,current_iq), g, dh)
     DO k = 1, dfftp%nnr
        vgg(k, is) = vgg(k, is) - dh(k)
     END DO
  END DO
  DO k = 1, dfftp%nnr
     dvxc(k,1) = 0.5d0*(vgg(k,1)+vgg(k,2))
     DO is = 2, 4
        term = SQRT(DOT_PRODUCT(drho(k,2:4), drho(k,2:4)))
        dvxc(k,is) = 0.5d0*segni(k)*( vgg(k,1)-vgg(k,2) &
        &                           + vsgga(k)*(drho(k,is)-term))
     END DO
  END DO
  !
  DEALLOCATE(dh)
  DEALLOCATE(h)
  DEALLOCATE(gdrho)
  DEALLOCATE(rhoout)
  DEALLOCATE(drhoout)
  DEALLOCATE(gdmag)
  DEALLOCATE(vgg)
  !
  RETURN
  !
END SUBROUTINE dgradcorr2
!
END MODULE sctk_dmuxc
