! mkl_stub.f90 -----------------------------------------------------------
! Stubs for the two Intel-MKL threading-control routines called by
! Callaway_2D_2V_DG.f90.  They are no-ops when we are not linking MKL.
!
! This file is *additive*: it does not modify any line of the reference
! solver, it only supplies symbols the reference expects at link time.
!-------------------------------------------------------------------------
SUBROUTINE MKL_SET_NUM_THREADS (N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER :: DUMMY
DUMMY = N
END SUBROUTINE MKL_SET_NUM_THREADS

SUBROUTINE MKL_SET_DYNAMIC (N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER :: DUMMY
DUMMY = N
END SUBROUTINE MKL_SET_DYNAMIC
