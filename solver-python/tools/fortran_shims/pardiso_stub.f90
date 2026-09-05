! pardiso_stub.f90 -------------------------------------------------------
! A drop-in replacement for the subset of the Intel MKL PARDISO API that
! Synthetic_Acceleration.f90 uses, implemented on top of reference LAPACK.
!
! Used only for the `--no-mkl` build flavour, so that the GSIS path of the
! Fortran reference can still be *run* on a machine without MKL.  The
! global trace system of the shipped case is N = N_FCS*3*NDOF_FC = 3840,
! which is small enough for a dense LU (118 MB, factorised exactly once).
!
! Supported phases:
!   12  analysis + numerical factorisation  (build dense matrix, DGETRF)
!   33  solve                               (DGETRS)
!   -1  release                             (deallocate)
! Anything else is a no-op.
!
! Additive only: no line of the reference solver is modified.
!-------------------------------------------------------------------------
MODULE PARDISO_STUB_STATE
IMPLICIT NONE
SAVE
INTEGER :: NSTORE = 0
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: LU
INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV
END MODULE PARDISO_STUB_STATE

SUBROUTINE PARDISOINIT (PT, MTYPE, IPARM)
IMPLICIT NONE
INTEGER*8, INTENT(OUT) :: PT(64)
INTEGER, INTENT(IN) :: MTYPE
INTEGER, INTENT(OUT) :: IPARM(64)
INTEGER :: I
DO I = 1, 64
   PT(I) = 0
   IPARM(I) = 0
END DO
IPARM(1) = 1
IPARM(3) = MTYPE   ! keep the argument used; value is irrelevant here
IPARM(3) = 1
END SUBROUTINE PARDISOINIT

SUBROUTINE PARDISO (PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, &
                    NRHS, IPARM, MSGLVL, B, X, ERROR)
USE PARDISO_STUB_STATE
IMPLICIT NONE
EXTERNAL :: DGETRF, DGETRS
INTEGER*8, INTENT(INOUT) :: PT(64)
INTEGER, INTENT(IN) :: MAXFCT, MNUM, MTYPE, PHASE, N, NRHS, MSGLVL
INTEGER, INTENT(IN) :: IA(*), JA(*), PERM(*)
INTEGER, INTENT(INOUT) :: IPARM(64)
DOUBLE PRECISION, INTENT(IN) :: A(*)
DOUBLE PRECISION, INTENT(INOUT) :: B(*)
DOUBLE PRECISION, INTENT(OUT) :: X(*)
INTEGER, INTENT(OUT) :: ERROR
INTEGER :: I, K, INFO, IDUM

ERROR = 0
IDUM = MAXFCT + MNUM + MTYPE + NRHS + MSGLVL + PERM(1) + IPARM(1)

IF (PHASE.EQ.12 .OR. PHASE.EQ.11 .OR. PHASE.EQ.22 .OR. PHASE.EQ.13) THEN
   IF (ALLOCATED(LU)) DEALLOCATE(LU)
   IF (ALLOCATED(IPIV)) DEALLOCATE(IPIV)
   ALLOCATE (LU(N,N), IPIV(N))
   LU = 0.d0
   DO I = 1, N
      DO K = IA(I), IA(I+1)-1
         LU(I,JA(K)) = LU(I,JA(K)) + A(K)
      END DO
   END DO
   CALL DGETRF (N, N, LU, N, IPIV, INFO)
   IF (INFO.NE.0) ERROR = -1
   NSTORE = N
   PT(1) = 1
END IF

IF (PHASE.EQ.33 .OR. PHASE.EQ.13) THEN
   IF (.NOT.ALLOCATED(LU) .OR. NSTORE.NE.N) THEN
      ERROR = -2
      RETURN
   END IF
   DO I = 1, N
      X(I) = B(I)
   END DO
   CALL DGETRS ('N', N, 1, LU, N, IPIV, X, N, INFO)
   IF (INFO.NE.0) ERROR = -3
END IF

IF (PHASE.EQ.-1 .OR. PHASE.EQ.0) THEN
   IF (ALLOCATED(LU)) DEALLOCATE(LU)
   IF (ALLOCATED(IPIV)) DEALLOCATE(IPIV)
   NSTORE = 0
   PT(1) = 0
END IF

END SUBROUTINE PARDISO
