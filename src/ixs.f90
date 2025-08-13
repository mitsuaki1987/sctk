MODULE ixs_mod
  IMPLICIT NONE
  !
  CONTAINS
SUBROUTINE form_factor(ispc,q,fq)
  !
  USE kinds, ONLY : dp
  USE constants, ONLY : pi
  !
  INTEGER,INTENT(IN) :: ispc
  REAL(dp),INTENT(IN) :: q ! [A^{-1}]
  REAL(dp),INTENT(OUT) :: fq
  !
  INTEGER :: ii
  REAL(dp) :: a(4,3), b(4,3), c(3) ! 1:Pd, 2:Sb, 3:Pt
  !Pd
  a(1,1) = 19.3319
  b(1,1) = 0.698655
  a(2,1) = 15.5017
  b(2,1) = 7.98929
  a(3,1) = 5.29537
  b(3,1) = 25.2052
  a(4,1) = 0.605844
  b(4,1) = 76.8986
  c(1) = 5.26593
  !Sb
  a(1,2) = 19.6418
  b(1,2) = 5.3034
  a(2,2) = 19.0455
  b(2,2) = 0.4607
  a(3,2) = 5.0371
  b(3,2) = 27.9074
  a(4,2) = 2.6827
  b(4,2) = 75.2825
  c(2) = 4.5909
  !Pt
  a(1,3) = 27.0059
  b(1,3) = 1.51293
  a(2,3) = 17.7639
  b(2,3) = 8.81174
  a(3,3) = 15.7131
  b(3,3) = 0.424593
  a(4,3) = 5.7837
  b(4,3) = 38.6103
  c(3) = 11.6883
  !
  fq = c(ispc)
  DO ii = 1, 4
    fq = fq + a(ii,ispc) * EXP(-b(ii,ispc)*(q/(4.0_dp*pi))**2)
  END DO
  !
END SUBROUTINE form_factor
!
FUNCTION bedist(beta,en) RESULT(dist)
  USE kinds, ONLY : dp
  REAL(dp),INTENT(IN) :: beta, en
  REAL(dp) :: dist
  dist = 0.5_dp *(1.0 / TANH(0.5_dp * beta * en) - 1.0_dp)
END FUNCTION bedist
END MODULE ixs_mod
!
PROGRAM ixs
  USE kinds, ONLY : dp
  USE constants, ONLY: pi, RY_TO_THZ, BOHR_RADIUS_ANGS, K_BOLTZMANN_RY
  USE ixs_mod, ONLY : form_factor, bedist
  IMPLICIT NONE
  !
  INTEGER :: nat, nmode, igv(3), ig1, ig2, ig3, itmp, iat, imode, iq, nq
  REAL(dp) :: bvec(3,3), qvec(3), qlen, rtmp, beta, alat, temperature, qlen2, dq(3)
  COMPLEX(dp) :: csq
  INTEGER,ALLOCATABLE :: spc(:)
  REAL(dp),ALLOCATABLE :: freq(:,:), mass(:), fq(:), dw(:), msd(:,:), q0vec(:,:), sq(:,:), qx(:)
  COMPLEX(dp),ALLOCATABLE :: eig(:,:,:,:), qeig(:)
  CHARACTER(100) :: filename
  !
  read(*,*) nat, temperature, nq
  temperature = temperature * K_BOLTZMANN_RY ! K -> Ry
  beta = 1.0_dp / temperature
  read(*,*) alat
  read(*,*) bvec(1:3,1)
  read(*,*) bvec(1:3,2)
  read(*,*) bvec(1:3,3)
  READ(*,*) filename
  bvec(1:3, 1:3) = bvec(1:3, 1:3) * 2.0_dp * pi / alat ! a.u.
  nmode = 3 * nat
  ALLOCATE(eig(3,nat,nmode,nq), freq(nmode,nq), mass(nat), spc(nat), fq(nat), qeig(nat), dw(nat), &
  &        msd(nat,nq), q0vec(3,nq), sq(nmode,nq), qx(nq))
  DO iat = 1, nat
    read(*,*) spc(iat), mass(iat)
  END DO
  !
  OPEN(100, file=TRIM(filename))
  !
  DO iq = 1, nq
    !
    READ(100, '(a)')
    READ(100, '(a)')
    READ(100, '(1x,6x,3f12.4)') q0vec(1:3,iq)
    READ(100, '(a)')
    q0vec(1:3,iq) = q0vec(1:3,iq) * 2.0_dp * pi / alat ! a.u.
    !
    do imode = 1, nmode
      !
      read(100,9010) itmp, freq(imode,iq), rtmp
      !
      do iat = 1,nat
        read(100,9020) eig(1:3,iat,imode,iq)
        eig(1:3,iat,imode,iq) = eig(1:3,iat,imode,iq) * SQRT(mass(iat))
      end do
      rtmp = REAL(SUM(eig(1:3,1:nat,imode,iq)*CONJG(eig(1:3,1:nat,imode,iq))), dp)
      eig(1:3,1:nat,imode,iq) = eig(1:3,1:nat,imode,iq) / SQRT(rtmp)
      !
    end do
    READ(100, '(a)')
    freq(1:nmode,iq) = freq(1:nmode,iq) / RY_TO_THZ ! THz -> Ry
    !
    ! Mean Square Displacement
    !
    DO iat = 1, nat
      msd(iat,iq) = 0.0_dp
      DO imode = 1, nmode
        msd(iat,iq) = msd(iat,iq) + REAL(DOT_PRODUCT(eig(1:3,iat,imode,iq), eig(1:3,iat,imode,iq)),dp) &
        &                   * (2.0_dp * bedist(beta,freq(imode,iq)) + 1.0_dp) / freq(imode,iq)
      END DO
      msd(iat,iq) = msd(iat,iq) / (3.0_dp * mass(iat))
    END DO
    !
  END DO !iq
  !
  qx(1) = 0.0_dp
  qlen2 = 9999.9
  DO iq = 2, nq
    dq(1:3) = q0vec(1:3,iq) - q0vec(1:3,iq-1)
    qlen = SQRT(DOT_PRODUCT(dq,dq))
    IF(qlen > qlen2*6) THEN
      qx(iq) = qx(iq-1)
    ELSE
      qx(iq) = qx(iq-1) + qlen
    END IF
    qlen2 = qlen
  END DO
  !
  CLOSE(100)
  !
  9010 format(5x, 6x, i5, 3x, f15.6, 8x, f15.6, 7x)
  9020 format(1x, 1x, 3(f10.6, 1x, f10.6, 3x), 1x)
  !
  DO ig1 = 0, 5
    DO ig2 = 0, 5
      DO ig3 = 0, 5
        !
        igv(1:3) = (/ig1, ig2, ig3/)
        !
        DO iq = 1, nq
          !
          qvec(1:3) = MATMUL(bvec(1:3,1:3), REAL(igv(1:3), dp)) + q0vec(1:3,iq)
          qlen = SQRT(DOT_PRODUCT(qvec, qvec))
          !
          DO iat = 1, nat
            CALL form_factor(spc(iat), qlen/BOHR_RADIUS_ANGS, fq(iat)) ! bohr^-1 -> ang^{-1}
            dw(iat) = EXP(-0.5_dp * qlen**2 *msd(iat,iq))
          END DO
          !
          DO imode = 1, nmode
            DO iat = 1, nat
              qeig(iat) = DOT_PRODUCT(qvec(1:3), eig(1:3,iat,imode,iq))
            END DO
            csq = SUM(fq(1:nat) * dw(1:nat) / SQRT(mass(1:nat)) * qeig(1:nat))
            sq(imode,iq) = (bedist(beta, freq(imode,iq)) + 1.0_dp) / freq(imode,iq) * REAL(csq*CONJG(csq),dp)
          END DO
          !
        END DO
        !
        WRITE(filename,'("ixs",i1,"-",i1,"-",i1,".dat")') ig1, ig2, ig3
        OPEN(200, file=TRIM(filename))
        DO imode = 1, nmode
          DO iq = 1, nq
            WRITE(200,*) qx(iq), freq(imode,iq)*RY_TO_THZ, sq(imode,iq)
          END DO
          WRITE(200,*)
        END DO
        CLOSE(200)
        !
      END DO
    END DO
  END DO
  !
END PROGRAM ixs
