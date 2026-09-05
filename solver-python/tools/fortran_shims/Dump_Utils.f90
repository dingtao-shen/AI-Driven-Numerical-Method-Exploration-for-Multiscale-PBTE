! Dump_Utils.f90 ---------------------------------------------------------
! Raw binary stage dumps for cross-validating the Python port against the
! Fortran reference (Proposal 1, Phase 1).
!
! Activated by the environment variable PYBTE_DUMP_DIR.  When it is unset
! every routine here is a cheap no-op, so the instrumented binary can also
! serve as the plain production reference.
!
! Layout produced in $PYBTE_DUMP_DIR:
!     manifest.txt      one line per array:  <name> <dtype> <ndim> <d1..dn>
!     <name>.bin        raw little-endian payload, Fortran (column-major)
!                       element order, exactly as the array sits in memory
!
! dtype is 'r8' (REAL*8) or 'i4' (default INTEGER).
!
! DUMP_R8/DUMP_I4 are deliberately *external* procedures with an implicit
! interface so that Fortran sequence association lets us hand them arrays
! of any rank without a per-rank generic or a RESHAPE temporary.
!
! Additive only: no line of the reference solver is modified by this file.
!-------------------------------------------------------------------------
MODULE DUMP_UTILS
IMPLICIT NONE
SAVE

LOGICAL :: DUMP_ON = .FALSE.
LOGICAL :: DUMP_VDF_ON = .TRUE.
CHARACTER (LEN=512) :: DUMP_DIR = ''
INTEGER, PARAMETER :: DUMP_MAX_ITERS = 64
INTEGER :: DUMP_ITERS(DUMP_MAX_ITERS) = 0
INTEGER :: DUMP_NITERS = 0
INTEGER, PARAMETER :: U_MANIFEST = 77
INTEGER, PARAMETER :: U_PAYLOAD  = 78

CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE DUMP_INIT ()
IMPLICIT NONE
INTEGER :: ios, ln, i, j, v
CHARACTER (LEN=256) :: SPEC

CALL GET_ENVIRONMENT_VARIABLE ('PYBTE_DUMP_DIR', DUMP_DIR, ln, ios)
IF (ios.NE.0 .OR. ln.EQ.0) THEN
   DUMP_ON = .FALSE.
   RETURN
END IF
DUMP_ON = .TRUE.

! iteration list, default 1,2,3,10
SPEC = ''
CALL GET_ENVIRONMENT_VARIABLE ('PYBTE_DUMP_ITERS', SPEC, ln, ios)
IF (ios.NE.0 .OR. ln.EQ.0) THEN
   DUMP_NITERS = 4
   DUMP_ITERS(1) = 1
   DUMP_ITERS(2) = 2
   DUMP_ITERS(3) = 3
   DUMP_ITERS(4) = 10
ELSE
   DUMP_NITERS = 0
   i = 1
   DO WHILE (i.LE.ln)
      j = i
      DO WHILE (j.LE.ln)
         IF (SPEC(j:j).EQ.',') EXIT
         j = j + 1
      END DO
      IF (j.GT.i) THEN
         READ (SPEC(i:j-1),*,IOSTAT=ios) v
         IF (ios.EQ.0 .AND. DUMP_NITERS.LT.DUMP_MAX_ITERS) THEN
            DUMP_NITERS = DUMP_NITERS + 1
            DUMP_ITERS(DUMP_NITERS) = v
         END IF
      END IF
      i = j + 1
   END DO
END IF

SPEC = ''
CALL GET_ENVIRONMENT_VARIABLE ('PYBTE_DUMP_NO_VDF', SPEC, ln, ios)
IF (ios.EQ.0 .AND. ln.GT.0) DUMP_VDF_ON = .FALSE.

! truncate the manifest
OPEN (U_MANIFEST, FILE=TRIM(DUMP_DIR)//'/manifest.txt', STATUS='REPLACE', ACTION='WRITE')
CLOSE (U_MANIFEST)

WRITE (*,*) '**** stage dumps enabled -> ', TRIM(DUMP_DIR), ' ****'

END SUBROUTINE DUMP_INIT

!-------------------------------------------------------------------------
LOGICAL FUNCTION DUMP_WANT_ITER (STEP)
IMPLICIT NONE
INTEGER, INTENT(IN) :: STEP
INTEGER :: i
DUMP_WANT_ITER = .FALSE.
IF (.NOT.DUMP_ON) RETURN
DO i = 1, DUMP_NITERS
   IF (DUMP_ITERS(i).EQ.STEP) THEN
      DUMP_WANT_ITER = .TRUE.
      RETURN
   END IF
END DO
END FUNCTION DUMP_WANT_ITER

!-------------------------------------------------------------------------
SUBROUTINE DUMP_MANIFEST_LINE (NAME, DTYPE, ND, DIMS)
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: NAME, DTYPE
INTEGER, INTENT(IN) :: ND
INTEGER, INTENT(IN) :: DIMS(ND)
INTEGER :: i
OPEN (U_MANIFEST, FILE=TRIM(DUMP_DIR)//'/manifest.txt', STATUS='UNKNOWN', &
      ACTION='WRITE', POSITION='APPEND')
WRITE (U_MANIFEST,'(A,1X,A,1X,I2)',ADVANCE='NO') TRIM(NAME), TRIM(DTYPE), ND
DO i = 1, ND
   WRITE (U_MANIFEST,'(1X,I12)',ADVANCE='NO') DIMS(i)
END DO
WRITE (U_MANIFEST,'(A)') ' '
CLOSE (U_MANIFEST)
END SUBROUTINE DUMP_MANIFEST_LINE

!-------------------------------------------------------------------------
SUBROUTINE DUMP_SCALAR_R8 (NAME, V)
IMPLICIT NONE
EXTERNAL :: DUMP_R8
CHARACTER (LEN=*), INTENT(IN) :: NAME
REAL (KIND=8), INTENT(IN) :: V
REAL (KIND=8) :: TMP(1)
INTEGER :: D(1)
TMP(1) = V
D(1) = 1
CALL DUMP_R8 (NAME, TMP, 1, 1, D)
END SUBROUTINE DUMP_SCALAR_R8

!-------------------------------------------------------------------------
SUBROUTINE DUMP_SCALAR_I4 (NAME, V)
IMPLICIT NONE
EXTERNAL :: DUMP_I4
CHARACTER (LEN=*), INTENT(IN) :: NAME
INTEGER, INTENT(IN) :: V
INTEGER :: TMP(1), D(1)
TMP(1) = V
D(1) = 1
CALL DUMP_I4 (NAME, TMP, 1, 1, D)
END SUBROUTINE DUMP_SCALAR_I4

!-------------------------------------------------------------------------
! "it000003_vdf_after_sweep"-style names
FUNCTION DUMP_ITNAME (STEP, TAG) RESULT (S)
IMPLICIT NONE
INTEGER, INTENT(IN) :: STEP
CHARACTER (LEN=*), INTENT(IN) :: TAG
CHARACTER (LEN=128) :: S
S = ''
WRITE (S,'(A2,I6.6,A1,A)') 'it', STEP, '_', TRIM(TAG)
END FUNCTION DUMP_ITNAME

END MODULE DUMP_UTILS

!=========================================================================
! External writers -- implicit interface on purpose (sequence association)
!=========================================================================
SUBROUTINE DUMP_R8 (NAME, X, NTOT, ND, DIMS)
USE DUMP_UTILS
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: NAME
INTEGER, INTENT(IN) :: NTOT, ND
INTEGER, INTENT(IN) :: DIMS(ND)
REAL (KIND=8), INTENT(IN) :: X(NTOT)
IF (.NOT.DUMP_ON) RETURN
OPEN (U_PAYLOAD, FILE=TRIM(DUMP_DIR)//'/'//TRIM(NAME)//'.bin', STATUS='REPLACE', &
      ACCESS='STREAM', FORM='UNFORMATTED', ACTION='WRITE')
WRITE (U_PAYLOAD) X
CLOSE (U_PAYLOAD)
CALL DUMP_MANIFEST_LINE (NAME, 'r8', ND, DIMS)
END SUBROUTINE DUMP_R8

SUBROUTINE DUMP_I4 (NAME, X, NTOT, ND, DIMS)
USE DUMP_UTILS
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: NAME
INTEGER, INTENT(IN) :: NTOT, ND
INTEGER, INTENT(IN) :: DIMS(ND)
INTEGER, INTENT(IN) :: X(NTOT)
IF (.NOT.DUMP_ON) RETURN
OPEN (U_PAYLOAD, FILE=TRIM(DUMP_DIR)//'/'//TRIM(NAME)//'.bin', STATUS='REPLACE', &
      ACCESS='STREAM', FORM='UNFORMATTED', ACTION='WRITE')
WRITE (U_PAYLOAD) X
CLOSE (U_PAYLOAD)
CALL DUMP_MANIFEST_LINE (NAME, 'i4', ND, DIMS)
END SUBROUTINE DUMP_I4
