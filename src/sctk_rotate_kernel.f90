!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_rotate_kernel
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Average matrix in grid
!
SUBROUTINE expand_g_v()
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp, ONLY : mp_bcast, mp_max
  USE symm_base, ONLY : nsym
  USE modes, ONLY : nmodes
  USE disp,  ONLY : nqs
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  USE sctk_val, ONLY : fbee, gg, gg0, lbee, nqbz, omg, omg0, Vc, Vc0, &
  &                    ngap1, kindx, lsf, nci
  !
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, ik, jk, ik2, ii, cnt, dsp, cntmax, ipe, fstk, lstk
  INTEGER,ALLOCATABLE :: qindx(:), nqsym(:,:), iks(:,:,:,:)
  REAL(dp),ALLOCATABLE :: gg1(:,:,:,:,:), Vc1(:,:,:,:,:), omg1(:,:)
  !
  CALL start_clock("expand_g_v")
  !
  CALL divide(world_comm, ngap1,ik,jk)
  fstk = MINVAL(kindx(ik:jk,1))
  lstk = MAXVAL(kindx(ik:jk,1))
  !
  ALLOCATE(nqsym(nqs, fstk:lstk), iks(2, nsym*2, nqs, fstk:lstk))
  !
  WRITE(stdout,'(9x,"Total RAM for Vc per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(lbee-fbee+1,dp)**2*REAL(nqbz,dp)*REAL(lstk-fstk+1,dp)*8.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for gg per process : ",e10.2," GB")') &
  &     REAL(nmodes,dp)*REAL(elph_nbnd_max-elph_nbnd_min+1,dp)**2 &
  &     *REAL(nqbz,dp)*REAL(lstk-fstk+1,dp)*8.0e-9_dp
  ALLOCATE(Vc(nci*lsf,fbee:lbee,nqbz,fbee:lbee,fstk:lstk), &
  &        gg(nmodes,elph_nbnd_min:elph_nbnd_max,nqbz,elph_nbnd_min:elph_nbnd_max,fstk:lstk), &
  &        omg(nmodes,nqbz,fstk:lstk))
  !
  CALL cnt_and_dsp(nqs,cnt,dsp)
  !
  cntmax = cnt
  CALL mp_max(cntmax, world_comm)
  !
  WRITE(stdout,'(9x,"Total RAM for Vc (temporary) per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(lbee-fbee+1,dp)**2*REAL(nqbz,dp)*REAL(cntmax,dp)*8.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for gg (temporary) per process : ",e10.2," GB")') &
  &     REAL(nmodes,dp)*REAL(elph_nbnd_max-elph_nbnd_min+1,dp)**2 &
  &     *REAL(nqbz,dp)*REAL(cntmax,dp)*8.0e-9_dp
  ALLOCATE(gg1(nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,nqbz,cntmax), &
  &        Vc1(nci*lsf,fbee:lbee,fbee:lbee,nqbz,cntmax), &
  &        omg1(nmodes,cntmax), qindx(cntmax))
  !
  ! #####  Symmetrize  #####
  !
  CALL k_kplusq_sym(fstk,lstk,nqsym,iks)
  !
  ! #####  Communicate  #####
  !
  Vc(1:nci*lsf,fbee:lbee,1:nqbz,fbee:lbee,fstk:lstk) = 0.0_dp
  gg( 1:nmodes,elph_nbnd_min:elph_nbnd_max,1:nqbz,elph_nbnd_min:elph_nbnd_max,fstk:lstk) = 0.0_dp
  omg(1:nmodes,          1:nqbz,          fstk:lstk) = 0.0_dp
  DO ipe = 0, nproc - 1
     !
     Vc1(1:nci*lsf,fbee:lbee,fbee:lbee,1:nqbz,1:cntmax) = 0.0_dp
     gg1( 1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz,1:cntmax) = 0.0_dp
     omg1(1:nmodes,                         1:cntmax) = 0.0_dp
     qindx(1:cntmax) = 1
     !
     IF(ipe == mpime) THEN
        !
        Vc1(1:nci*lsf,fbee:lbee,fbee:lbee,1:nqbz,1:cnt) = Vc0(1:nci*lsf,fbee:lbee,fbee:lbee,1:nqbz,1:cnt)
        gg1(     1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz,1:cnt) &
        & = gg0( 1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,1:nqbz,1:cnt)
        omg1(1:nmodes,                         1:cnt) = omg0(1:nmodes,                         1:cnt)
        !        
        DO iq = 1, cnt
           qindx(iq) = dsp + iq
        END DO
        !
     END IF
     !
     CALL mp_bcast(Vc1,   ipe, world_comm )
     CALL mp_bcast(gg1,   ipe, world_comm )
     CALL mp_bcast(omg1,  ipe, world_comm )
     CALL mp_bcast(qindx, ipe, world_comm )
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(nqsym,iks,fstk,lstk,cntmax,qindx,lsf,nmodes,nci, &
     !$OMP &        elph_nbnd_min, elph_nbnd_max, Vc1,Vc,omg1,omg,gg1,gg,fbee,lbee) &
     !$OMP & PRIVATE(ik,jk,iq,ik2,ii)
     !
     !$OMP DO
     DO ik = fstk, lstk
        DO iq = 1, cntmax
           !
           DO ii = 1, nqsym(qindx(iq),ik)
              !
              jk  = iks(1, ii, qindx(iq), ik)
              ik2 = iks(2, ii, qindx(iq), ik)
              !
              Vc(1:nci*lsf,fbee:lbee,jk,fbee:lbee,ik) = Vc(1:nci*lsf,fbee:lbee,jk,fbee:lbee,ik) &
              &                   + Vc1(1:nci*lsf,fbee:lbee,fbee:lbee,ik2,iq)
              !
              gg(    1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,elph_nbnd_min:elph_nbnd_max,ik) &
              & = gg(1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,elph_nbnd_min:elph_nbnd_max,ik) &
              &        + gg1(1:nmodes,elph_nbnd_min:elph_nbnd_max,elph_nbnd_min:elph_nbnd_max,ik2,iq)
              !
              omg(1:nmodes,jk,ik) = omg(1:nmodes,jk,ik) &
              &        + omg1(1:nmodes,iq)
              !
           END DO ! ii
           !
        END DO ! ik2
     END DO ! iq
     !$OMP END DO
     !$OMP END PARALLEL
     !
  END DO ! ipe
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,lsf,nmodes,nqs,nqbz,Vc,gg,omg,iks,nsym, &
  !$OMP &        elph_nbnd_min,elph_nbnd_max,fbee,lbee,nci) &
  !$OMP & PRIVATE(ik,jk)
  !
  !$OMP DO
  DO ik = fstk, lstk
     DO jk = 1, nqbz
        Vc(1:nci*lsf,fbee:lbee,jk,fbee:lbee,ik) = Vc(1:nci*lsf,fbee:lbee,jk,fbee:lbee,ik) &
        &                         / REAL(count(iks(1,1:nsym*2,1:nqs,ik) == jk), dp)
        gg(    1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,elph_nbnd_min:elph_nbnd_max,ik) &
        & = gg(1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,elph_nbnd_min:elph_nbnd_max,ik) &
        &                                   / REAL(count(iks(1,1:nsym*2,1:nqs,ik) == jk), dp)
        omg(1:nmodes,jk,ik) = omg(1:nmodes,jk,ik) &
        &       / REAL(count(iks(1,1:nsym*2,1:nqs,ik) == jk), dp)
     END DO ! jk
  END DO ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  DEALLOCATE(gg0, Vc0, omg0, gg1, Vc1, omg1, nqsym, iks)
  !
  CALL stop_clock("expand_g_v")
  !
END SUBROUTINE expand_g_v
!
! Symmetrize
!
SUBROUTINE k_kplusq_sym(fstk,lstk,nqsym,iks)
  !
  USE kinds, ONLY : DP
  USE symm_base, ONLY : nsym, s, time_reversal
  USE cell_base, ONLY : at
  USE disp,  ONLY : nq1, nq2, nq3, nqs, x_q
  USE sctk_val, ONLY : nqbz
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: fstk, lstk
  INTEGER,INTENT(OUT) :: nqsym(nqs,fstk:lstk)
  INTEGER,INTENT(OUT) :: iks(2,nsym*2,nqs,fstk:lstk)
  !
  INTEGER :: iq, ik, isym, ikv1(3), ik2, jk2, tsign(2), ntsign, itsign
  REAL(dp) :: kv0(3), kv1(3), rqv(3)
  !
  tsign(1:2) = (/1, -1/)
  if(time_reversal) THEN
     ntsign = 2
  ELSE
     ntsign = 1
  END IF
  !
  nqsym(           1:nqs,fstk:lstk) = 0
  iks(1:2,1:nsym*2,1:nqs,fstk:lstk) = 0
  !
  DO iq = 1, nqs
     !
     rqv(1:3) = MATMUL(x_q(1:3,iq), at(1:3,1:3))
     !
     DO ik = 1, nqbz
        !
        ikv1(1) = (ik - 1) / (nq2*nq3)
        ikv1(2) = (ik - 1 - ikv1(1)*nq2*nq3) / nq3
        ikv1(3) =  ik - 1 - ikv1(1)*nq2*nq3 - ikv1(2)*nq3
        kv0(1:3) = REAL(ikv1(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
        !
        DO isym = 1, nsym
           !
           DO itsign = 1, ntsign
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), kv0(1:3)) &
              &        * REAL((/nq1,nq2,nq3/), dp)
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              ik2 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              IF(ik2 < fstk .OR. lstk < ik2) CYCLE
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), kv0(1:3) + rqv(1:3))
              kv1(1:3) = kv1(1:3) * REAL((/nq1,nq2,nq3/), dp) - 0.5_dp
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              jk2 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              nqsym(iq,ik2) = nqsym(iq,ik2) + 1
              iks(1:2,nqsym(iq,ik2),iq,ik2) = (/jk2, ik/)
              !
           END DO ! itsign
           !
        END DO ! isym = 1, nsym
        !
     END DO ! ik
     !
  END DO ! isym
  !
END SUBROUTINE k_kplusq_sym
!
! Average matrix in grid
!
SUBROUTINE expand_g_v_f()
  !
  USE kinds, ONLY : DP
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp, ONLY : mp_bcast, mp_max
  USE io_global, ONLY : stdout
  USE modes, ONLY : nmodes
  USE klist, ONLY : nks
  USE start_k, ONLY : nk1, nk2, nk3
  USE disp,  ONLY : nqs
  USE fermisurfer_common, ONLY : b_low, b_high
  USE el_phon, ONLY : elph_nbnd_min, elph_nbnd_max
  !
  USE sctk_val, ONLY : fbee, gg0, ggF, lbee, nqbz, omg0, omgF, Vc0, VcF, lsf, nci
  USE sctk_cnt_dsp, ONLY : cnt_and_dsp
  !
  INTEGER :: iq, ik, jk, ik2, ii, cnt, dsp, cntmax, ipe, nindmax, &
  &          nqsym(nqs), nind(nqbz,nqbz), ikv(3), nks0, nks1
  REAL(dp) :: kv(3)
  INTEGER,ALLOCATABLE :: qindx(:), iks(:,:,:), ind(:,:,:,:)
  REAL(dp),ALLOCATABLE :: gg1(:,:,:,:,:), Vc1(:,:,:,:,:), omg1(:,:), wght(:,:), wght0(:,:)
  !
  CALL start_clock("expand_g_v_f")
  !
  CALL divide(world_comm, nks,nks0,nks1)
  CALL cnt_and_dsp(nqs, cnt, dsp)
  !
  ! #####  Expand with Symmetries  #####
  !
  CALL k_kplusq_sym_f(nindmax,nind,ind)
  !
  ! #####  Communicate  #####
  !
  cntmax = cnt
  CALL mp_max(cntmax, world_comm)
  !
  WRITE(stdout,'(9x,"Total RAM for Vc per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(lbee-fbee+1,dp)*REAL(b_high-b_low+1,dp)&
  &     *REAL(nqbz,dp)*REAL(nks1-nks0+1,dp)*8.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for gg per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(elph_nbnd_max-elph_nbnd_min+1,dp)*REAL(b_high-b_low+1,dp)&
  &     *REAL(nqbz,dp)*REAL(nks1-nks0+1,dp)*8.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for Vc (temporary) per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(lbee-fbee+1,dp)*REAL(b_high-b_low+1,dp)&
  &     *REAL(nqbz,dp)*REAL(cntmax,dp)*8.0e-9_dp
  WRITE(stdout,'(9x,"Total RAM for gg (temporary) per process : ",e10.2," GB")') &
  &     REAL(nci*lsf,dp)*REAL(elph_nbnd_max-elph_nbnd_min+1,dp)*REAL(b_high-b_low+1,dp)&
  &     *REAL(nqbz,dp)*REAL(cntmax,dp)*8.0e-9_dp
  ALLOCATE(Vc1(nci*lsf,         fbee:lbee,         b_low:b_high,nqbz,cntmax), &
  &        gg1( nmodes,elph_nbnd_min:elph_nbnd_max,b_low:b_high,nqbz,cntmax), &
  &       omg1( nmodes,                                              cntmax), &
  &      qindx(                                                      cntmax+1))
  !
  ALLOCATE(VcF(nci*lsf,fbee:lbee,                  nqbz,b_low:b_high,nks0:nks1), &
  &        ggF( nmodes,elph_nbnd_min:elph_nbnd_max,nqbz,b_low:b_high,nks0:nks1), &
  &       omgF( nmodes,                            nqbz,             nks0:nks1), &
  &      wght0(                                    nqbz,             nks0:nks1)  )
  !
  VcF(1:nci*lsf,            fbee:lbee,      1:nqbz,b_low:b_high,nks0:nks1) = 0.0_dp
  ggF( 1:nmodes,elph_nbnd_min:elph_nbnd_max,1:nqbz,b_low:b_high,nks0:nks1) = 0.0_dp
  omgF(1:nmodes,                            1:nqbz,             nks0:nks1) = 0.0_dp
  wght0(                                    1:nqbz,             nks0:nks1) = 0.0_dp
  DO ipe = 0, nproc - 1
     !
     Vc1(1:nci*lsf,         fbee:lbee,         b_low:b_high,1:nqbz,1:cntmax) = 0.0_dp
     gg1( 1:nmodes,elph_nbnd_min:elph_nbnd_max,b_low:b_high,1:nqbz,1:cntmax) = 0.0_dp
     omg1(1:nmodes,                                                1:cntmax) = 0.0_dp
     qindx(                                                        1:cntmax) = 1
     !
     IF(ipe == mpime) THEN
        !
        Vc1(    1:nci*lsf,        fbee:lbee,         b_low:b_high,1:nqbz,1:cnt) &
        & = Vc0(1:nci*lsf,        fbee:lbee,         b_low:b_high,1:nqbz,1:cnt)
        gg1(    1:nmodes,elph_nbnd_min:elph_nbnd_max,b_low:b_high,1:nqbz,1:cnt) &
        & = gg0(1:nmodes,elph_nbnd_min:elph_nbnd_max,b_low:b_high,1:nqbz,1:cnt)
        omg1(   1:nmodes,                                                1:cnt) &
        & = omg0(1:nmodes,                                               1:cnt)
        !
        qindx(cntmax + 1) = cnt
        DO iq = 1, cnt
           qindx(iq) = dsp + iq
        END DO
        !
     END IF
     !
     CALL mp_bcast(Vc1,   ipe, world_comm)
     CALL mp_bcast(gg1,   ipe, world_comm)
     CALL mp_bcast(omg1,  ipe, world_comm)
     CALL mp_bcast(qindx, ipe, world_comm)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(nks0,nks1,cntmax,qindx,lsf,fbee,lbee,nmodes, &
     !$OMP &        Vc1,VcF,omg1,omgF,gg1,ggF,nind,nindmax,ind,wght0,nk1,nk2,nk3, &
     !$OMP &        b_low,b_high,nci,elph_nbnd_min,elph_nbnd_max) &
     !$OMP & PRIVATE(ik,jk,iq,ik2,ii,iks,wght,nqsym,ikv,kv)
     !
     !$OMP DO
     DO ik = nks0, nks1
        !
        ! #####  Interpol & Symmetrize  #####
        !
        ikv(1) = (ik - 1) / (nk2*nk3)
        ikv(2) = (ik - 1 - ikv(1)*nk2*nk3) / nk3
        ikv(3) =  ik - 1 - ikv(1)*nk2*nk3 - ikv(2)*nk3
        kv(1:3) = REAL(ikv(1:3), dp) / REAL((/nk1,nk2,nk3/), dp)
        CALL interpol_g_v(kv,nind,nindmax,ind,nqsym,iks,wght)
        !
        DO iq = 1, qindx(cntmax + 1)
           !
           DO ii = 1, nqsym(qindx(iq))
              !
              jk  = iks(1, ii, qindx(iq))
              ik2 = iks(2, ii, qindx(iq))
              !
              wght0(jk,ik) = wght0(jk,ik) + wght(ii, qindx(iq))
              !
              VcF(    1:nci*lsf,fbee:lbee,jk,b_low:b_high,ik) &
              & = VcF(1:nci*lsf,fbee:lbee,jk,b_low:b_high,ik) &
              & + Vc1(1:nci*lsf,fbee:lbee,   b_low:b_high,ik2,iq) * wght(ii, qindx(iq))
              !
              ggF(    1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,b_low:b_high,ik) &
              & = ggF(1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,b_low:b_high,ik) &
              & + gg1(1:nmodes,elph_nbnd_min:elph_nbnd_max,   b_low:b_high,ik2,iq) * wght(ii, qindx(iq))
              !
              omgF(1:nmodes,jk,ik) = omgF(1:nmodes,jk,ik) &
              &                    + omg1(1:nmodes,iq) * wght(ii, qindx(iq))
              !
           END DO ! ii
           !
        END DO ! iq
        !
        DEALLOCATE(iks,wght)
        !
     END DO ! ik
     !$OMP END DO
     !$OMP END PARALLEL
     !
  END DO ! ipe
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nks0,nks1,lsf,nmodes,fbee,lbee,nqbz,VcF,ggF,omgF,wght0, &
  !$OMP &        elph_nbnd_min,elph_nbnd_max,b_low,b_high,nci) &
  !$OMP & PRIVATE(ik,jk)
  !
  !$OMP DO
  DO ik = nks0, nks1
     DO jk = 1, nqbz
        VcF(     1:nci*lsf,        fbee:lbee,         jk,b_low:b_high,ik) &
        & = VcF( 1:nci*lsf,        fbee:lbee,         jk,b_low:b_high,ik) / wght0(jk,ik)
        ggF(     1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,b_low:b_high,ik) &
        & = ggF( 1:nmodes,elph_nbnd_min:elph_nbnd_max,jk,b_low:b_high,ik) / wght0(jk,ik)
        omgF(    1:nmodes,                            jk,             ik) &
        & = omgF(1:nmodes,                            jk,             ik) / wght0(jk,ik)
     END DO ! jk
  END DO ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  DEALLOCATE(gg0, Vc0, omg0)
  DEALLOCATE(gg1, Vc1, omg1, ind, wght0)
  !
  CALL stop_clock("expand_g_v_f")
  !
END SUBROUTINE expand_g_v_f
!
! Average
!
SUBROUTINE k_kplusq_sym_f(nindmax,nind,ind)
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE symm_base, ONLY : nsym, s, time_reversal
  USE cell_base, ONLY : at
  USE disp,  ONLY : nq1, nq2, nq3, nqs, x_q
  USE sctk_val, ONLY : nqbz
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(OUT) :: nind(nqbz,nqbz), nindmax
  INTEGER,INTENT(OUT),ALLOCATABLE :: ind(:,:,:,:)
  !
  INTEGER :: ik, iq, isym, jk1, ik1, ikv1(3), ntsign, itsign, tsign(2)
  REAL(dp) :: kv0(3), kv1(3), rqv(3)
  !
  tsign(1:2) = (/1, -1/)
  if(time_reversal) THEN
     ntsign = 2
  ELSE
     ntsign = 1
  END IF
  !
  ! Query of memory size
  !
  nind(1:nqbz,1:nqbz) = 0
  !
  DO ik = 1, nqbz
     !
     ikv1(1) = (ik - 1) / (nq2*nq3)
     ikv1(2) = (ik - 1 - ikv1(1)*nq2*nq3) / nq3
     ikv1(3) =  ik - 1 - ikv1(1)*nq2*nq3 - ikv1(2)*nq3
     kv0(1:3) = REAL(ikv1(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
     !   
     DO iq = 1, nqs
        !
        rqv(1:3) = MATMUL(x_q(1:3,iq), at(1:3,1:3))
        !
        DO isym = 1, nsym
           !
           DO itsign = 1, ntsign
              !
              ! rotate k
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), &
              &                 kv0(1:3)             ) * REAL((/nq1,nq2,nq3/), dp)
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              ik1 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              ! Rotate k + q
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp),  &
              &                 kv0(1:3) + rqv(1:3)) * REAL((/nq1,nq2,nq3/), dp)
              kv1(1:3) = kv1(1:3) - 0.5_dp
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              jk1 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              nind(jk1,ik1) = nind(jk1,ik1) + 1
              !
           END DO ! itsign
           !
        END DO ! isym
        !
     END DO ! iq
     !
  END DO ! ik
  !
  nindmax = MAXVAL(nind(1:nqbz,1:nqbz))
  ALLOCATE(ind(2,nindmax,nqbz,nqbz))
  !
  ind(1:2,1:nindmax,1:nqbz,1:nqbz) = 0
  nind(1:nqbz,1:nqbz) = 0
  !
  WRITE(stdout,'(7x,"nindmax : ",i0)') nindmax
  !
  DO ik = 1, nqbz
     !
     ikv1(1) = (ik - 1) / (nq2*nq3)
     ikv1(2) = (ik - 1 - ikv1(1)*nq2*nq3) / nq3
     ikv1(3) =  ik - 1 - ikv1(1)*nq2*nq3 - ikv1(2)*nq3
     kv0(1:3) = REAL(ikv1(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
     !
     DO iq = 1, nqs
        !
        rqv(1:3) = MATMUL(x_q(1:3,iq), at(1:3,1:3))
        !
        DO isym = 1, nsym
           !
           DO itsign = 1, ntsign
              !
              ! rotate k
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), &
              &                 kv0(1:3)             ) * REAL((/nq1,nq2,nq3/), dp)
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              ik1 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              ! Rotate k + q
              !
              kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp),  &
              &                 kv0(1:3) + rqv(1:3)) * REAL((/nq1,nq2,nq3/), dp)
              kv1(1:3) = kv1(1:3) - 0.5_dp
              ikv1(1:3) = NINT(kv1(1:3))
              !
              IF(ANY(ABS(kv1(1:3) - REAL(ikv1(1:3), dp)) > 1e-5_dp)) CYCLE
              !
              ikv1(1:3) = MODULO(ikv1(1:3), (/nq1,nq2,nq3/))
              jk1 = 1 + ikv1(3) + ikv1(2)*nq3 + ikv1(1)*nq2*nq3
              !
              nind(jk1,ik1) = nind(jk1,ik1) + 1
              ind(1:2,nind(jk1,ik1),jk1,ik1) = (/ik, iq/)
              !
           END DO ! itsign
           !
        END DO ! isym
        !
     END DO ! iq
     !
  END DO ! ik
  !
END SUBROUTINE k_kplusq_sym_f
!
! Symmetrize
!
SUBROUTINE interpol_g_v(kv0,nind,nindmax,ind,nqsym,iks,wght)
  !
  USE kinds, ONLY : DP
  USE symm_base, ONLY : nsym, s, time_reversal
  USE disp,  ONLY : nq1, nq2, nq3, nqs
  USE sctk_val, ONLY : nqbz
  !
  USE sctk_tetra, ONLY : interpol_indx
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nind(nqbz,nqbz), nindmax, ind(2,nindmax,nqbz,nqbz)
  REAL(dp),INTENT(IN) :: kv0(3)
  INTEGER,INTENT(OUT) :: nqsym(nqs)
  INTEGER,INTENT(OUT),ALLOCATABLE :: iks(:,:,:)
  REAL(dp),INTENT(OUT),ALLOCATABLE :: wght(:,:)
  !
  INTEGER :: ii, jj, jk, isym, ik2, jk1, iq, jkv1(3), ikintp(8), nqsmax, &
  &          tsign(2), itsign, ntsign
  REAL(dp) :: wintp(8), kv1(3)
  !
  tsign(1:2) = (/1, -1/)
  if(time_reversal) THEN
     ntsign = 2
  ELSE
     ntsign = 1
  END IF
  !
  nqsym(1:nqs) = 0
  !
  ! Query of memory size
  !
  DO isym = 1, nsym
     !
     DO itsign = 1, ntsign
        !
        kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), kv0(1:3))
        CALL interpol_indx((/nq1,nq2,nq3/),kv1,ikintp,wintp)
        !
        DO jk = 1, nqbz
           !
           jkv1(1) = (jk - 1) / (nq2*nq3)
           jkv1(2) = (jk - 1 - jkv1(1)*nq2*nq3) / nq3
           jkv1(3) =  jk - 1 - jkv1(1)*nq2*nq3 - jkv1(2)*nq3
           kv1(1:3) = REAL(jkv1(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
           !
           kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), &
           &                 kv1(1:3) + 0.5_dp / REAL((/nq1,nq2,nq3/), dp))
           kv1(1:3) = kv1(1:3) * REAL((/nq1,nq2,nq3/), dp) - 0.5_dp
           jkv1(1:3) = NINT(kv1(1:3))
           !
           IF(ANY(ABS(kv1(1:3) - REAL(jkv1(1:3), dp)) > 1e-5_dp)) CYCLE
           !
           jkv1(1:3) = MODULO(jkv1(1:3), (/nq1,nq2,nq3/))
           jk1 = 1 + jkv1(3) + jkv1(2)*nq3 + jkv1(1)*nq2*nq3
           !
           DO ii = 1, 8
              DO jj = 1, nind(jk1,ikintp(ii))
                 !
                 iq = ind(2,jj,jk1,ikintp(ii))
                 !
                 nqsym(iq) = nqsym(iq) + 1
                 !
              END DO ! jj
           END DO ! ii
           !
        END DO ! jk
        !
     END DO ! itsign
     !
  END DO ! isym
  !
  nqsmax = MAXVAL(nqsym(1:nqs))
  !
  ALLOCATE(wght(nqsmax,nqs), iks(2,nqsmax,nqs))
  !
  iks(1:2,1:nqsmax,1:nqs) = 0
  wght(   1:nqsmax,1:nqs) = 0.0_dp
  nqsym(             1:nqs) = 0
  !
  DO isym = 1, nsym
     !
     DO itsign = 1, ntsign
        !
        kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), kv0(1:3))
        CALL interpol_indx((/nq1,nq2,nq3/),kv1,ikintp,wintp)
        !
        DO jk = 1, nqbz
           !
           jkv1(1) = (jk - 1) / (nq2*nq3)
           jkv1(2) = (jk - 1 - jkv1(1)*nq2*nq3) / nq3
           jkv1(3) =  jk - 1 - jkv1(1)*nq2*nq3 - jkv1(2)*nq3
           kv1(1:3) = REAL(jkv1(1:3), dp) / REAL((/nq1,nq2,nq3/), dp)
           !
           kv1(1:3) = MATMUL(REAL(s(1:3,1:3,isym)*tsign(itsign), dp), &
           &                 kv1(1:3) + 0.5_dp / REAL((/nq1,nq2,nq3/), dp))
           kv1(1:3) = kv1(1:3) * REAL((/nq1,nq2,nq3/), dp) - 0.5_dp
           jkv1(1:3) = NINT(kv1(1:3))
           !
           IF(ANY(ABS(kv1(1:3) - REAL(jkv1(1:3), dp)) > 1e-5_dp)) CYCLE
           !
           jkv1(1:3) = MODULO(jkv1(1:3), (/nq1,nq2,nq3/))
           jk1 = 1 + jkv1(3) + jkv1(2)*nq3 + jkv1(1)*nq2*nq3
           !
           DO ii = 1, 8
              DO jj = 1, nind(jk1,ikintp(ii))
                 !
                 ik2 = ind(1,jj,jk1,ikintp(ii))
                 iq  = ind(2,jj,jk1,ikintp(ii))
                 !
                 nqsym(iq) = nqsym(iq) + 1
                 iks( 1:2, nqsym(iq), iq) = (/jk, ik2/)
                 wght(     nqsym(iq), iq) = wintp(ii)
                 !
              END DO ! jj
           END DO ! ii
           !
        END DO ! jk
        !
     END DO ! itsign = 1, ntsign
     !
  END DO ! isym
  !
END SUBROUTINE interpol_g_v
  !
END MODULE sctk_rotate_kernel
