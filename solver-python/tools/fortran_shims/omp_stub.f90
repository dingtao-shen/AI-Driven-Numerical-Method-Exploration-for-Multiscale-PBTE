! omp_stub.f90 -----------------------------------------------------------
! The reference driver calls OMP_SET_NUM_THREADS(24) unconditionally.
! We build the validation reference *serially* (no -fopenmp) so that the
! residual history is bit-reproducible, which means the symbol has to come
! from somewhere.  This is a no-op.
!
! Additive only: no line of the reference solver is modified.
!-------------------------------------------------------------------------
SUBROUTINE OMP_SET_NUM_THREADS (N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER :: DUMMY
DUMMY = N
END SUBROUTINE OMP_SET_NUM_THREADS
