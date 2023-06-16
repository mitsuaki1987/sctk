!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_z
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Calc Z_{n k}
!
SUBROUTINE make_Z()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime, world_comm
  USE mp, ONLY : mp_sum
  USE modes, ONLY : nmodes
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : beta, bindx, dk, emin, gg, lz_coulomb, freq_min, &
  &                    kindx, ngap, ngap1, ngap2, omg, xi, Z, lsf, Vc, nci
  !
  USE sctk_kernel_weight, ONLY : Zweight, calc_Zsf
  !
  IMPLICIT NONE
  !
  INTEGER :: igap, ik, ib, jgap, jk, jb, imode, ngap10, ngap11
  REAL(dp) :: x, xp, om, zave(2), ave(2), Z0
  !
  CALL start_clock("make_Z")
  !
  IF(.NOT. ALLOCATED(Z)) ALLOCATE(Z(ngap,2))
  Z(1:ngap,1:2) = 0.0_dp
  !
  CALL divide(world_comm, ngap1,ngap10,ngap11)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ngap10,ngap11,ngap2,nmodes,kindx,bindx,beta,lz_coulomb, &
  !$OMP &        elph_nbnd_min,elph_nbnd_max,xi,dk,Z,gg,omg,lsf,nci,Vc,freq_min) &
  !$OMP & PRIVATE(igap,ik,ib,jgap,jk,jb,imode,x,xp,om,Z0)
  !
  !$OMP DO REDUCTION(+: Z)
  DO igap = ngap10, ngap11
     !
     x = xi(igap,1)
     ik = kindx(igap,1)
     ib = bindx(igap,1)
     !
     DO jgap = 1, ngap2
        !
        xp = xi(jgap,2)
        jk = kindx(jgap,2)
        jb = bindx(jgap,2)
        !
        ! Spin-fluctuation renormalization
        !
        IF(lsf==2 .OR. lz_coulomb) THEN
           !
           IF(lsf == 2) THEN
              Z0 =   calc_Zsf(ABS(x)+ABS(xp), Vc(nci+1:nci*2,jb,jk,ib,ik))
           ELSE
              Z0 = - calc_Zsf(ABS(x)+ABS(xp), Vc(1:nci,jb,jk,ib,ik))
           END IF
           !
           Z(igap,1) = Z(igap,1) + dk(jgap,2) * Z0
           Z(jgap,2) = Z(jgap,2) + dk(igap,1) * Z0
           !
        END IF
        !      
        ! Electron-phonon renormalization
        !
        IF(     elph_nbnd_min <= ib .AND. ib <= elph_nbnd_max &
        & .AND. elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
           !
           DO imode = 1, nmodes
              !
              om = omg(imode,jk,ik)
              IF(om < freq_min) CYCLE
              !
              Z(igap,1) = Z(igap,1) - dk(jgap,2) * gg(imode,jb,jk,ib,ik) &
              &                     * Zweight(beta, ABS(x), ABS(xp), om)
              Z(jgap,2) = Z(jgap,2) - dk(igap,1) * gg(imode,jb,jk,ib,ik) &
              &                     * Zweight(beta, ABS(xp), ABS(x), om)
              !
           END DO ! imode
           !
        END IF ! (elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max)
        !
     END DO ! jgap
     !
  END DO ! igap
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  CALL mp_sum( Z, world_comm )
  !
  zave(1) = SUM(Z(1:ngap,1) * dk(1:ngap,1), ABS(xi(1:ngap,1)) < emin)
  zave(2) = SUM(Z(1:ngap,2) * dk(1:ngap,2), ABS(xi(1:ngap,2)) < emin)
  ave(1) = SUM(dk(1:ngap,1), ABS(xi(1:ngap,1)) < emin)
  ave(2) = SUM(dk(1:ngap,2), ABS(xi(1:ngap,2)) < emin)
  !
  zave(1:2) = zave(1:2) / ave(1:2)
  !
  IF(mpime == 0)  WRITE(*,'(7x,"Averaged Z_{FS} : ",2(e12.5,2x))') zave(1:2)
  !
  CALL stop_clock("make_Z")
  !
END SUBROUTINE make_Z
!
! Calc Z_{n k}
!
SUBROUTINE make_Z_qpdos()
  !
  USE kinds, ONLY : DP
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE fermisurfer_common, ONLY : b_low, b_high
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : beta, bindx, dk, ggf, kindx, ngap2, nx, freq_min, &
  &                    omgf, xi, xi0, ZF, lsf, Vcf, nci, lz_coulomb
  USE sctk_kernel_weight, ONLY : Zweight, calc_Zsf
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, imode, ix, nks0, nks1
  REAL(dp) :: x, xp, om
  !
  CALL start_clock("make_Z_qpdos")
  !
  ALLOCATE(ZF(nx,b_low:b_high,nks))
  ZF(1:nx,b_low:b_high,1:nks) = 0.0_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nx,nks0,nks1,b_low,b_high,nmodes,ngap2,xi,dk,kindx,bindx,lsf,nci, &
  !$OMP &        ZF,ggf,omgf,xi0,elph_nbnd_min,elph_nbnd_max,beta,VcF,lz_coulomb,freq_min) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,ix,imode,x,xp,om)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        DO ix = 1, nx
           !
           x = ABS(xi0(ix))
           !
           DO jgap = 1, ngap2
              !
              xp = xi(jgap,2)
              jk = kindx(jgap,2)
              jb = bindx(jgap,2)
              !
              ! Spin-fluctuation renormalization
              !
              IF(lsf==2) THEN
                 ZF(ix,ib,ik) = ZF(ix,ib,ik) &
                 &            + dk(jgap,2) * calc_Zsf(ABS(x)+ABS(xp), Vcf(nci+1:nci*2,jb,jk,ib,ik))
              ELSE IF(lz_coulomb) THEN
                 ZF(ix,ib,ik) = ZF(ix,ib,ik) &
                 &            - dk(jgap,2) * calc_Zsf(ABS(x)+ABS(xp), Vcf(1:nci,jb,jk,ib,ik))
              END IF
              !
              ! Electron-phonon renormalization
              !
              IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
                 !
                 DO imode = 1, nmodes
                    !
                    om = omgf(imode,jk,ik)
                    IF(om < freq_min) CYCLE
                    !
                    ZF(ix,ib,ik) = ZF(ix,ib,ik) &
                    &        - dk(jgap,2) * ggf(imode,jb,jk,ib,ik) &
                    &        * Zweight(beta, ABS(x), ABS(xp), om)
                    !
                 END DO ! imode
                 !
              END IF ! (elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max)
              !
           END DO ! jgap
           !
        END DO
        !
     END DO ! ib
     !
  END DO ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  CALL mp_sum( ZF, world_comm )
  !
  CALL stop_clock("make_Z_qpdos")
  !
END SUBROUTINE make_Z_qpdos
!
! Calc Z_{n k}
!
SUBROUTINE make_Z_f()
  !
  USE kinds, ONLY : DP
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE mp_world, ONLY : world_comm
  USE mp, ONLY : mp_sum
  USE fermisurfer_common, ONLY : b_low, b_high
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : beta, bindx, dk, ggf, kindx, ngap2, freq_min, &
  &                    omgf, xi, ZF, lsf, Vcf, nci, lz_coulomb
  USE sctk_kernel_weight, ONLY : Zweight, calc_Zsf
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib, jgap, jk, jb, imode, nks0, nks1
  REAL(dp) :: xp, om
  !
  CALL start_clock("make_Z_f")
  !
  ALLOCATE(ZF(1,b_low:b_high,nks))
  ZF(1,b_low:b_high,1:nks) = 0.0_dp
  !
  CALL divide(world_comm, nks,nks0,nks1)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,b_low,b_high,nmodes,ngap2,xi,dk,kindx,bindx,ZF,nci, &
  !$OMP &        ggf,omgf,elph_nbnd_min,elph_nbnd_max,beta,lsf,Vcf,lz_coulomb,freq_min) &
  !$OMP & PRIVATE(ik,ib,jgap,jk,jb,imode,xp,om)
  !
  !$OMP DO
  DO ik = nks0, nks1
     !
     DO ib = b_low, b_high
        !
        DO jgap = 1, ngap2
           !
           xp = xi(jgap,2)
           jk = kindx(jgap,2)
           jb = bindx(jgap,2)
           !
           ! Spin-fluctuation renormalization
           !
           IF(lsf==2) THEN
              ZF(1,ib,ik) = ZF(1,ib,ik) &
              &           + dk(jgap,2) * calc_Zsf(ABS(xp), Vcf(nci+1:nci*2,jb,jk,ib,ik))
           ELSE IF(lz_coulomb) THEN
              ZF(1,ib,ik) = ZF(1,ib,ik) &
              &           - dk(jgap,2) * calc_Zsf(ABS(xp), Vcf(    1:nci,jb,jk,ib,ik))
           END IF
           !
           ! Electron-phonon renormalization
           !
           IF(elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max) THEN
              !
              DO imode = 1, nmodes
                 !
                 om = omgf(imode,jk,ik)
                 IF(om < freq_min) CYCLE
                 !
                 ZF(1,ib,ik) = ZF(1,ib,ik) &
                 &         - dk(jgap,2) * ggf(imode,jb,jk,ib,ik) &
                 &         * Zweight(beta, 0.0_dp, ABS(xp), om)
                 !
              END DO ! imode
              !
           END IF ! (elph_nbnd_min <= jb .AND. jb <= elph_nbnd_max)
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
  CALL mp_sum( ZF, world_comm )
  !
  CALL stop_clock("make_Z_f")
  !
END SUBROUTINE make_Z_f
!
END MODULE sctk_z
