!
! Copyright (C) 2017 Mitsuaki Kawamura
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE sctk_val
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & lsf,           & !< =2 for spin-fluctuation, =0 for non-sf
  & fbee,          & !< First band for Vc
  & laddxc,        & !< Switch for XC
  & lbee,          & !< Last band for Vc
  & nmf,           & !< Number of Matsubara frequencies
  & nci, &
  & ne,            & !<
  & nqbz,          & !< nq1 * nq2 * nq3
  & nf(3),         & !< # grid for FFT
  & nftot,         & !< Total # of G
  & ngap,          & !< Total # of Delta, Xi, ...
  & ngap1,         & !< Total # of Delta, Xi, ... for grid 1 (w/o k-shift)
  & ngap2,         & !< Total # of Delta, Xi, ... for grid 2 (with k-shift)
  & ngv,           & !< # of G-vector in ecutwfc
  & ngv0,          & !< # of G-vector in ecutwfc
  & ngv1,          & !< # of G-vector in ecutwfc
  & nkpe,          & !< # of k par PE
  & nx               !< = 2 * ne - 1
  !
  REAL(8),SAVE :: &
  & beta,       & !< inversed temperature [Ry]
  & emax,       & !< Max energy for qpdos
  & emin,       & !< Minimum energy scale [Ry]
  & Zemin,      & !< Minimum energy scale [Ry]
  & xic           !< Cut off Xi
  !
  LOGICAL,SAVE :: &
  & zero_kelvin
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & bindx(:,:),   & !< (ngap,2) band index for gap equation
  & gindx(:),     & !< (nftot) G-vector in ecutwfc
  & igv(:,:,:,:), & !< (3,npwmax,nk,2). G points
  & kindx(:,:),   & !< (ngap,2) k point for gap equation
  & npw(:,:)        !< (nk,2). # of PWs
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & delta(:,:),     & !< (ngap,2) Kohn-Sham gap functions [Ry]
  & dk(:,:),        & !< (ngap,2) Weight of k
  & dltF(:,:,:),    & !< (nx,nbf,nk) Kohn-Sham gap functions [Ry] at Fermi surface
  & dosk(:,:,:),    & !< (nx,nbnd,nk)
  & dx0(:),         & !< (nx) weight for energy
  & e0(:),          & !< (ne) Energy grid for qpdos
  & effint(:,:),    & !< (ngap1,ngap2) Effective interaction    
  & Fvel(:,:,:),    & !< (3,b_low:b_high, nk) The Fermi velocity
  & gg(:,:,:,:,:),  & !< (nm,nbnd,nks,nbnd,nqbz*) El-Ph matrix element [Ry]
  & gg0(:,:,:,:,:), & !< (nm,nbnd,nbnd,nqbz,nqs*) El-Ph matrix element [Ry]
  & ggF(:,:,:,:,:), & !< (nm,nbnd,nqbz,nbf,nks) El-Ph matrix element [Ry]
  & gq2(:),         & !< (nftot) |G+q|^2
  & Kel(:,:,:,:,:),  & !< (0:nci,nbnd,nbnd,nqbz,2). Coulomb Kernel
  & Ksf(:,:,:,:),   & !< (2,nbnd,nbnd,nqbz). Spin-fluctuation Kernel
  & mf(:),          & !< (nmf) for SCDFT, (nci) for Kel. Matsubara frequencies
  & omg(:,:,:),     & !< (nmodes,nqbz,nqbz*) Phonon frequencies [Ry]
  & omg0(:,:),      & !< (nmodes,nqs) Phonon frequencies [Ry]
  & omgF(:,:,:),    & !< (nmodes,nqbz,nks) Phonon frequencies [Ry]
  & sdos(:),        & !< (ne) Superconducting QPDSO
  & Vc(:,:,:,:,:),  & !< (nci,nb,nk,nb,nk*) Screened Coulomb matrix element [Ry]
  & Vc0(:,:,:,:,:), & !< (nci,nb,nb,nk,nk0) Screened Coulomb matrix element [Ry]
  & VcF(:,:,:,:,:), & !< (nci,nb,nk,nbf,nks) Screened Coulomb matrix element [Ry]
  & wmf(:,:),       & !< (nmf) Weights for Matsubara frequency
  & xi(:,:),        & !< (ngap,2) Kohn-Sham energy [Ry]
  & xi0(:),         & !< (nx) energy scale [Ry]
  & Z(:,:),         & !< (ngap,2) Renormalization factor for grid 1      
  & ZF(:,:,:)         !< (nx,nbf,nks) Renormalization factor ar Fermi surface
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & wfc(:,:,:,:,:), & !< (npwmax,npol,nbnd,nks,2) wfc(G)
  & wfc1(:,:,:,:),  & !< (nfto,npol,nbnd,nks) wfc(ib,ik)
  & wfc1q(:,:,:,:), & !< (nftot,npol,nbnd,nks) wfc(ib,ik)
  & wfc2(:,:,:,:),  & !< (nftot,npol,nbnd,nks) wfc(jb,jk)
  & wfc2q(:,:,:,:), & !< (nftot,npol,nbnd,nks) wfc(jb,jk)
  & wscr(:,:,:,:)     !< (nftot,nftot,0:nci,2*npol) Screened interaction
  !
END MODULE sctk_val
