MODULE ACCELERATION
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE INTEGRATIONS
USE BOUNDARY_CONDITION
USE VELOCITY_GRID
USE VELOCITY_DISTRIBUTION_FUNCTION
IMPLICIT NONE

REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: UQ !DIMENSION(7*NDOF_TRI,N_TRIS), macro properties and derivatives of velocity and temperature
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: U_TRACE !DIMENSION(3*NDOF_FC,N_FCS), traces of macro properties on interior faces
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:), SAVE :: ST !DIMENSION(3), stabilization factors
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:,:), SAVE :: inv_AA_SOL !DIMENSION(7*NDOF_TRI, 7*NDOF_TRI, N_TRIS) coefficient matrix of UQ in local problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:,:,:), SAVE :: AA_TRACE !DIMENSION(7*NDOF_TRI, 3*NDOF_FC, 3, N_TRIS) coefficient matrix of traces in local problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: AA_SRC !DIMENSION(7*NDOF_TRI, N_TRIS) source term

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: KKA ! DIMENSION(N_FCS*NDOF_FC*3, N_FCS*NDOF_FC*3), global matrix in global problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: AA_TMP, D1A_TMP, D2A_TMP
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: BA_SOL !DIMENSION(3*NDOF_FC,7*NDOF_TRI,3,N_TRIS), coefficient matrix for the field variables in global problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: KKAcomp, C1A_TMP, C2A_TMP
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: IKKA, JKKA, JOBA

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: FFA, HHA !R.H.S of global problem, DIMENSION(N_FCS*5*NDOF_FC)

!***********************************************************************************************
!pardiso
INTEGER, SAVE :: MAXFCT, MTYPE, PHASE, ERROR, MSGLVL
INTEGER*8, SAVE :: PT(64)
INTEGER, SAVE :: IPARM(64)
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: PERM

CONTAINS

SUBROUTINE Init_Acceleration_New ()
IMPLICIT NONE
EXTERNAL :: DGETRF, DGETRI
INTEGER :: I, M, L, IL, LWORK, INFO, TRID1, TRID2, LFCID1, LFCID2, GFCID1, GFCID2, GFCID,BCID, PAIR, NNZ
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: AA, BB, CC, DD, OO, PP, WW
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: AA_SOL
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:) :: WORK_A, KKAcomp_tmp
INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV_A, JKKA_tmp
REAL (KIND=DBL) :: nx, ny

ALLOCATE(UQ(7*NDOF_TRI,N_TRIS),U_TRACE(3*NDOF_FC,N_FCS),ST(3))
UQ = 0.d0
U_TRACE = 0.d0
ST(1) = 1.d0
ST(2) = 1.d0!/Hmin!/DELTA_REF
ST(3) = 1.d0!/Hmin!/DELTA_REF

ALLOCATE (inv_AA_SOL(7*NDOF_TRI,7*NDOF_TRI,N_TRIS))
inv_AA_SOL = 0.d0

ALLOCATE(AA(NDOF_TRI,NDOF_TRI),BB(NDOF_TRI,NDOF_TRI),CC(NDOF_TRI,NDOF_TRI))
ALLOCATE(DD(NDOF_TRI,NDOF_TRI),OO(NDOF_TRI,NDOF_TRI),PP(NDOF_TRI,NDOF_TRI))
ALLOCATE(WORK_A(7*NDOF_TRI),IPIV_A(7*NDOF_TRI),AA_SOL(7*NDOF_TRI,7*NDOF_TRI))

LWORK = 7*NDOF_TRI

!$OMP PARALLEL DO PRIVATE(I,M,L,IL,INFO,AA,BB,CC,DD,PP,OO,AA_SOL,WORK_A,IPIV_A)
DO I = 1, N_TRIS
    AA = INT_NODFUNC_TRI_TRI_X(I,:,:)
    BB = INT_NODFUNC_TRI_TRI_Y(I,:,:)
    DO M = 1, NDOF_TRI
        DO L = 1, NDOF_TRI
            CC(M,L) = AA(L,M)
            DD(M,L) = BB(L,M)
        END DO
    END DO
    PP = 0.d0
    DO IL = 1, 3
        PP = PP + INT_NODFUNC_TRI_FC(I,IL,:,:)
    END DO
    OO = INT_NODFUNC_TRI_TRI(:,:,I)

    AA_SOL = 0.d0
    !T
    AA_SOL(1:NDOF_TRI,            1:  NDOF_TRI) = ST(1)*PP
    AA_SOL(1:NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = CC
    AA_SOL(1:NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = DD

    !qx
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI,            1:  NDOF_TRI) = -AA*Cv/3.d0
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = ST(2)*PP+OO/TAU_R
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 3*NDOF_TRI+1:4*NDOF_TRI) = -4.d0/3.d0*CC
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 4*NDOF_TRI+1:5*NDOF_TRI) = -DD
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 5*NDOF_TRI+1:6*NDOF_TRI) = -DD
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 6*NDOF_TRI+1:7*NDOF_TRI) =  2.d0/3.d0*CC

    !qy
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI,            1:  NDOF_TRI) = -BB*Cv/3.d0
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = ST(3)*PP+OO/TAU_R
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 3*NDOF_TRI+1:4*NDOF_TRI) =  2.d0/3.d0*DD
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 4*NDOF_TRI+1:5*NDOF_TRI) = -CC
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 5*NDOF_TRI+1:6*NDOF_TRI) = -CC
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 6*NDOF_TRI+1:7*NDOF_TRI) = -4.d0/3.d0*DD


    !Lxx
    AA_SOL(3*NDOF_TRI+1:4*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = AA*TAU_C/5.d0
    AA_SOL(3*NDOF_TRI+1:4*NDOF_TRI, 3*NDOF_TRI+1:4*NDOF_TRI) = OO

    !Lxy
    AA_SOL(4*NDOF_TRI+1:5*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = BB*TAU_C/5.d0
    AA_SOL(4*NDOF_TRI+1:5*NDOF_TRI, 4*NDOF_TRI+1:5*NDOF_TRI) = OO

    !Lyx
    AA_SOL(5*NDOF_TRI+1:6*NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = AA*TAU_C/5.d0
    AA_SOL(5*NDOF_TRI+1:6*NDOF_TRI, 5*NDOF_TRI+1:6*NDOF_TRI) = OO

    !Lyy
    AA_SOL(6*NDOF_TRI+1:7*NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = BB*TAU_C/5.d0
    AA_SOL(6*NDOF_TRI+1:7*NDOF_TRI, 6*NDOF_TRI+1:7*NDOF_TRI) = OO

    WORK_A = 0.d0
    IPIV_A = 0
    CALL DGETRF (LWORK,LWORK,AA_SOL,LWORK,IPIV_A,INFO)
    IF (INFO.NE.0) WRITE (*,*) 'Warning LU factorization error: ',INFO
    CALL DGETRI (LWORK,AA_SOL,LWORK,IPIV_A,WORK_A,LWORK,INFO)
    IF (INFO.NE.0) WRITE (*,*) 'Warning Inverse Matrix error: ',INFO
    inv_AA_SOL(:,:,I) = AA_SOL

    !DO M=1,NDOF_TRI
    !DO L=1,NDOF_TRI
    !    WRITE(*,*)PP(M,L),OO(M,L)
    !END DO
    !END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(AA,BB,CC,DD,PP,OO,AA_SOL,WORK_A,IPIV_A)

!*******************************************************************************
ALLOCATE(AA_TRACE(7*NDOF_TRI,3*NDOF_FC,3,N_TRIS))
ALLOCATE(WW(NDOF_TRI,NDOF_FC))
AA_TRACE = 0.d0

!$OMP PARALLEL DO PRIVATE(I,IL,WW,nx,ny)
DO I = 1, N_TRIS
    DO IL = 1, 3
        WW = INT_NODFUNC_TRI_FC_FC(:,:,I,IL)
        nx = TRIANGLES_INF(I,2*IL)
        ny = TRIANGLES_INF(I,2*IL+1)

        !T
        AA_TRACE(1:NDOF_TRI,           1:  NDOF_FC, IL, I) = ST(1)*WW

        !qx
        AA_TRACE(NDOF_TRI+1:2*NDOF_TRI,           1:  NDOF_FC, IL, I) = -nx*WW*Cv/3.d0
        AA_TRACE(NDOF_TRI+1:2*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = ST(2)*WW

        !qy
        AA_TRACE(2*NDOF_TRI+1:3*NDOF_TRI,           1:  NDOF_FC, IL, I) = -ny*WW*Cv/3.d0
        AA_TRACE(2*NDOF_TRI+1:3*NDOF_TRI, 2*NDOF_FC+1:3*NDOF_FC, IL, I) = ST(3)*WW

    
        !Lxx
        AA_TRACE(3*NDOF_TRI+1:4*NDOF_TRI, NDOF_FC+1:2*NDOF_FC, IL, I) = nx*WW*TAU_C/5.d0

        !Lxy
        AA_TRACE(4*NDOF_TRI+1:5*NDOF_TRI, NDOF_FC+1:2*NDOF_FC, IL, I) = ny*WW*TAU_C/5.d0

        !Lyx
        AA_TRACE(5*NDOF_TRI+1:6*NDOF_TRI, 2*NDOF_FC+1:3*NDOF_FC, IL, I) = nx*WW*TAU_C/5.d0

        !Lyy
        AA_TRACE(6*NDOF_TRI+1:7*NDOF_TRI, 2*NDOF_FC+1:3*NDOF_FC, IL, I) = ny*WW*TAU_C/5.d0

        AA_TRACE(:,:,IL,I) = MATMUL(inv_AA_SOL(:,:,I),AA_TRACE(:,:,IL,I))
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(WW)

!**********************************************************************************
ALLOCATE(AA_SRC(7*NDOF_TRI, N_TRIS))
AA_SRC = 0.d0

!***********************************************************************************
ALLOCATE(KKA(NDOF_FC*3,N_FCS*NDOF_FC*3))
ALLOCATE(BA_SOL(3*NDOF_FC,7*NDOF_TRI,3,N_TRIS))
ALLOCATE(KKAcomp_tmp(5*N_FCS*3*NDOF_FC*3*NDOF_FC))
ALLOCATE(IKKA(N_FCS*3*NDOF_FC+1))
ALLOCATE(JKKA_tmp(5*N_FCS*3*NDOF_FC*3*NDOF_FC))
ALLOCATE(FFA(N_FCS*3*NDOF_FC),HHA(N_FCS*3*NDOF_FC))

ALLOCATE(BB(NDOF_FC,NDOF_TRI))

KKA = 0.d0
BA_SOL = 0.d0
IKKA = 0
JKKA_tmp = 0
KKAcomp_tmp = 0.d0
FFA = 0.d0
HHA = 0.d0

!*****************************************************************************
!$OMP PARALLEL DO PRIVATE(I,IL,M,L,BB,nx,ny,GFCID,BCID)
DO I = 1, N_TRIS
    DO IL = 1, 3
        DO M = 1, NDOF_FC
            DO L = 1, NDOF_TRI
                BB(M,L) = INT_NODFUNC_TRI_FC_FC(L,M,I,IL)
            END DO
        END DO
        nx = TRIANGLES_INF(I,2*IL)
        ny = TRIANGLES_INF(I,2*IL+1)
        GFCID = TRIANGLES_TAG(I,3+IL)
        BCID = FACES_TAG(GFCID,3)

        IF (BCID.EQ.0) THEN
            !interior &
            BA_SOL(          1:  NDOF_FC,           1:  NDOF_TRI, IL, I) = 0.5d0*BB
            BA_SOL(          1:  NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = nx/2.d0/ST(1)*BB
            BA_SOL(          1:  NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = ny/2.d0/ST(1)*BB

            BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = 0.5d0*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = -2.d0*nx/3.d0/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -ny/2.d0/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -ny/2.d0/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) =  nx/3.d0/ST(2)*BB

            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = 0.5d0*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) =  ny/3.d0/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -nx/2.d0/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -nx/2.d0/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) = -2.d0*ny/3.d0/ST(3)*BB
        END IF
        IF (BCID.GT.0) THEN
            IF (BC_TYP(BCID).EQ.3) THEN
                !Periodic BC
                BA_SOL(          1:  NDOF_FC,           1:  NDOF_TRI, IL, I) = 0.5d0*BB
                BA_SOL(          1:  NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = nx/2.d0/ST(1)*BB
                BA_SOL(          1:  NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = ny/2.d0/ST(1)*BB

                BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = 0.5d0*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = -2.d0*nx/3.d0/ST(2)*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -ny/2.d0/ST(2)*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -ny/2.d0/ST(2)*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) =  nx/3.d0/ST(2)*BB

                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = 0.5d0*BB
                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) =  ny/3.d0/ST(3)*BB
                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -nx/2.d0/ST(3)*BB
                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -nx/2.d0/ST(3)*BB
                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) = -2.d0*ny/3.d0/ST(3)*BB
            END IF
            IF (BC_TYP(BCID).EQ.1) THEN
            !thermarlising wall
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = -4.d0*nx/3.d0/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -ny/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -ny/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) =  2.d0*nx/3.d0/ST(2)*BB

            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) =  2.d0*ny/3.d0/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = -nx/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = -nx/ST(3)*BB
            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,6*NDOF_TRI+1:7*NDOF_TRI, IL, I) = -4.d0*ny/3.d0/ST(3)*BB
            END IF
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(BB)

ALLOCATE(D1A_TMP(3*NDOF_FC,3*NDOF_FC),D2A_TMP(3*NDOF_FC,3*NDOF_FC))
ALLOCATE(AA(NDOF_FC,NDOF_FC),AA_TMP(3*NDOF_FC,3*NDOF_FC))
ALLOCATE(C1A_TMP(3*NDOF_FC),C2A_TMP(3*NDOF_FC))

NNZ = 0
DO I = 1, N_FCS
    KKA = 0.d0

    AA = INT_NODFUNC_FC_FC(:,:,I)
    AA_TMP(          1:  NDOF_FC,          1:  NDOF_FC) = AA
    AA_TMP(  NDOF_FC+1:2*NDOF_FC,  NDOF_FC+1:2*NDOF_FC) = AA
    AA_TMP(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_FC+1:3*NDOF_FC) = AA

    BCID = FACES_TAG(I,3)

    IF (BCID.EQ.0) THEN
        !interior faces
        TRID1 = FACES_TAG(I,5)
        TRID2 = FACES_TAG(I,6)
        LFCID1 = FACES_TAG(I,7)
        LFCID2 = FACES_TAG(I,8)

        KKA(1:3*NDOF_FC,1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = AA_TMP
        DO IL = 1,3
            !interior
            GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
            GFCID2 = TRIANGLES_TAG(TRID2,3+IL)

            D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
            KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) = &
                KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) - D1A_TMP

            D2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_TRACE(:,:,IL,TRID2))
            KKA(1:3*NDOF_FC, 1+(GFCID2-1)*3*NDOF_FC:GFCID2*3*NDOF_FC) = &
                KKA(1:3*NDOF_FC, 1+(GFCID2-1)*3*NDOF_FC:GFCID2*3*NDOF_FC) - D2A_TMP
        END DO
    ELSE
        !boundary
        IF (BC_TYP(BCID).EQ.3) THEN
            PAIR = FACES_TAG(I,4)
            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            IF (FACES_TAG(PAIR,5).EQ.0) THEN
                TRID2 = FACES_TAG(PAIR,6)
                LFCID2 = FACES_TAG(PAIR,8)
            ELSE
                TRID2 = FACES_TAG(PAIR,5)
                LFCID2 = FACES_TAG(PAIR,7)
            END IF

            KKA(1:3*NDOF_FC,1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = AA_TMP
            DO IL = 1,3
                !interior
                GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
                GFCID2 = TRIANGLES_TAG(TRID2,3+IL)

                D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
                KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) = &
                    KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) - D1A_TMP

                D2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_TRACE(:,:,IL,TRID2))
                KKA(1:3*NDOF_FC, 1+(GFCID2-1)*3*NDOF_FC:GFCID2*3*NDOF_FC) = &
                    KKA(1:3*NDOF_FC, 1+(GFCID2-1)*3*NDOF_FC:GFCID2*3*NDOF_FC) - D2A_TMP
            END DO

        ELSE

            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF

            KKA(1:3*NDOF_FC,1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = AA_TMP
            DO IL = 1,3
                !interior
                GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
                D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
                KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) = &
                    KKA(1:3*NDOF_FC, 1+(GFCID1-1)*3*NDOF_FC:GFCID1*3*NDOF_FC) - D1A_TMP
            END DO
        END IF
    END IF

    !compression
    DO M = 1, 3*NDOF_FC
        IKKA(M+(I-1)*3*NDOF_FC) = NNZ !row: M+(I-1)*3*NDOF_FC
        DO L = 1, N_FCS*3*NDOF_FC
            IF (ABS(KKA(M,L)).GT.1.d-30) THEN
                NNZ = NNZ + 1
                IF (NNZ.GT.(5*N_FCS*3*NDOF_FC*3*NDOF_FC)) THEN
                    WRITE(*,*)'Range Exceed'
                    STOP
                END IF
                KKAcomp_tmp(NNZ) = KKA(M,L)
                JKKA_tmp(NNZ) = L
            END IF
        END DO
    END DO
END DO

IKKA(N_FCS*3*NDOF_FC+1) = NNZ
IKKA = IKKA + 1
ALLOCATE(KKAcomp(NNZ),JKKA(NNZ))
KKAcomp = KKAcomp_tmp(1:NNZ)
JKKA = JKKA_tmp(1:NNZ)


DEALLOCATE(KKA,AA,AA_TMP,D1A_TMP,D2A_TMP,KKAcomp_tmp,JKKA_tmp)

END SUBROUTINE

SUBROUTINE Calculate_SRC_ACC_HoTfromDVM ()
IMPLICIT NONE
REAL (KIND=DBL) :: PIxx, PIxy, PIyy, VDF_tmp
INTEGER :: I, M, L, J1, J2

AA_SRC = 0.d0

!$OMP PARALLEL DO PRIVATE(I,M,PIxx,PIxy,PIyy,VDF_tmp,J1,J2,L)
DO I = 1, N_TRIS
    DO M = 1, NDOF_TRI
        PIxx = 0.d0
        PIxy = 0.d0
        PIyy = 0.d0

        !*****************************************************************
        DO J2 = 1, NAZIM
        DO J1 = 1, NPOLE
            VDF_tmp = 0.d0
            DO L = 1, NDOF_TRI
                !partial x
                VDF_tmp = VDF_tmp + VDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_X(I,L,M)
            END DO

            PIxx = PIxx + CX(J1,J2)*(5.d0*CX(J1,J2)*CX(J1,J2)-3.d0)*VDF_tmp*DOMEGA(J1,J2)
            PIxy = PIxy + CY(J1,J2)*(5.d0*CX(J1,J2)*CX(J1,J2)-1.d0)*VDF_tmp*DOMEGA(J1,J2)
            PIyy = PIyy + CX(J1,J2)*(5.d0*CY(J1,J2)*CY(J1,J2)-1.d0)*VDF_tmp*DOMEGA(J1,J2)

            VDF_tmp = 0.d0
            DO L = 1, NDOF_TRI
                !partial y
                VDF_tmp = VDF_tmp + VDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_Y(I,L,M)
            END DO

            PIxx = PIxx + CY(J1,J2)*(5.d0*CX(J1,J2)*CX(J1,J2)-1.d0)*VDF_tmp*DOMEGA(J1,J2)
            PIxy = PIxy + CX(J1,J2)*(5.d0*CY(J1,J2)*CY(J1,J2)-1.d0)*VDF_tmp*DOMEGA(J1,J2)
            PIyy = PIyy + CY(J1,J2)*(5.d0*CY(J1,J2)*CY(J1,J2)-3.d0)*VDF_tmp*DOMEGA(J1,J2)

        END DO
        END DO

        AA_SRC(M+3*NDOF_TRI,I) = (PIxx + 0.5d0*PIyy)*TAU_C/5.d0
        AA_SRC(M+4*NDOF_TRI,I) = 0.5d0*PIxy*TAU_C/5.d0
        AA_SRC(M+5*NDOF_TRI,I) = 0.5d0*PIxy*TAU_C/5.d0
        AA_SRC(M+6*NDOF_TRI,I) = (0.5d0*PIxx + PIyy)*TAU_C/5.d0
    END DO
    AA_SRC(:,I) = MATMUL(inv_AA_SOL(:,:,I),AA_SRC(:,I))
END DO
!$OMP END PARALLEL DO

END SUBROUTINE


SUBROUTINE Global_Problem_Solver_ACC ()
IMPLICIT NONE
EXTERNAL :: PARDISO
INTEGER :: I,TRID1,TRID2,LFCID1,LFCID2,M,L,J1,J2,BCID,PAIR
REAL (KIND=DBL) :: nx, ny
REAL (KIND=DBL) :: Tbc, qxbc, qybc, VDF_tmp


U_TRACE = 0.d0
FFA = 0.d0

!$OMP PARALLEL DO PRIVATE(I,BCID,TRID1,TRID2,LFCID1,LFCID2,C1A_TMP,C2A_TMP,nx,ny,M,qxbc,qybc,Tbc, &
!$OMP L,VDF_tmp,J1,J2,PAIR)
DO I = 1, N_FCS 
    BCID = FACES_TAG(I,3)
    IF (BCID.EQ.0) THEN
        !interior
        TRID1 = FACES_TAG(I,5)
        TRID2 = FACES_TAG(I,6)
        LFCID1 = FACES_TAG(I,7)
        LFCID2 = FACES_TAG(I,8)

        C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
        C2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_SRC(:,TRID2))

        FFA(1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = C1A_TMP + C2A_TMP
    ELSE
        !boundary
        IF (BC_TYP(BCID).EQ.3) THEN
            !periodic
            PAIR = FACES_TAG(I,4)
            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            IF (FACES_TAG(PAIR,5).EQ.0) THEN
                TRID2 = FACES_TAG(PAIR,6)
                LFCID2 = FACES_TAG(PAIR,8)
            ELSE
                TRID2 = FACES_TAG(PAIR,5)
                LFCID2 = FACES_TAG(PAIR,7)
            END IF
            C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
            C2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_SRC(:,TRID2))

            FFA(1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = C1A_TMP + C2A_TMP
        ELSE

            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            nx = TRIANGLES_INF(TRID1,2*LFCID1)
            ny = TRIANGLES_INF(TRID1,2*LFCID1+1)

            C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
            FFA(1+(I-1)*3*NDOF_FC:I*3*NDOF_FC) = C1A_TMP

            IF (BC_TYP(BCID).EQ.2) THEN
                !adiabatic WALL
                DO M = 1, NDOF_FC
                    qxbc = 0.d0
                    qybc = 0.d0
                    Tbc = 0.d0

                    DO J2 = 1, NAZIM
                    DO J1 = 1, NPOLE
                        VDF_tmp = 0.d0
                        DO L = 1, NDOF_TRI
                            VDF_tmp = VDF_tmp + VDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
                        END DO
                        Tbc  = Tbc  +           VDF_tmp*DOMEGA(J1,J2)
                        qxbc = qxbc + CX(J1,J2)*VDF_tmp*DOMEGA(J1,J2)
                        qybc = qybc + CY(J1,J2)*VDF_tmp*DOMEGA(J1,J2)
                    END DO
                    END DO
                    FFA(M+          (I-1)*3*NDOF_FC) = Tbc/Cv
                    FFA(M+  NDOF_FC+(I-1)*3*NDOF_FC) = qxbc !qxbc*ny*ny - qybc*nx*ny
                    FFA(M+2*NDOF_FC+(I-1)*3*NDOF_FC) = qybc !-qxbc*nx*ny + qybc*nx*nx
                END DO
            END IF

            IF (BC_TYP(BCID).EQ.1) THEN
                !thermarlising wall
                DO M = 1, NDOF_FC
                    Tbc = 0.d0

                    DO J2 = 1, NAZIM
                    DO J1 = 1, NPOLE
                        VDF_tmp = 0.d0
                        DO L = 1, NDOF_TRI
                            VDF_tmp = VDF_tmp + VDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
                        END DO
                        Tbc  = Tbc  +           VDF_tmp*DOMEGA(J1,J2)
                    END DO
                    END DO
                    FFA(M+          (I-1)*3*NDOF_FC) = Tbc/Cv
                END DO
                !DO M = 1, NDOF_FC
                !    FFA(M+(I-1)*3*NDOF_FC) = BC_TEMP(BCID)*INT_NODFUNC_FC(M,I)
                !END DO
            END IF
        END IF
    END IF
END DO
!$OMP END PARALLEL DO

!******************************************************************************************************
M = N_FCS*3*NDOF_FC
PHASE = 33
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,M,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)
IF (ERROR.NE.0) WRITE(*,*) 'PARDISO solution error', ERROR

!$OMP PARALLEL DO PRIVATE(I)
DO I = 1, N_FCS
    U_TRACE(:,I) = HHA(1+(I-1)*3*NDOF_FC:I*3*NDOF_FC)
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Local_Problem_Solver_ACC ()
IMPLICIT NONE
INTEGER :: I, IL

!$OMP PARALLEL DO PRIVATE (I,IL)
DO I = 1, N_TRIS
    UQ(:,I) = AA_SRC(:,I)
    DO IL = 1,3
        UQ(:,I) = UQ(:,I) + MATMUL(AA_TRACE(:,:,IL,I),U_TRACE(:,TRIANGLES_TAG(I,3+IL)))
    END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Init_PARDISO ()
IMPLICIT NONE
EXTERNAL :: pardisoinit, PARDISO
INTEGER :: N

WRITE(*,*) 'Initialize PARDISO ...'
MAXFCT = 1
MTYPE = 11
MSGLVL = 0

N = N_FCS*3*NDOF_FC
ALLOCATE(PERM(N))
PERM = 0

CALL pardisoinit(PT,MTYPE,IPARM)

IPARM(27) = 1
PHASE = 12
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,N,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)
IF (ERROR.NE.0) WRITE(*,*) 'PARDISO analysis & factorization error', ERROR
IF (ERROR.EQ.0) WRITE(*,*) 'PARDISO LU factorization successful'

END SUBROUTINE

SUBROUTINE Correct_VDF_Calculate_Macro_Properties ()
IMPLICIT NONE
INTEGER :: I, M, J1, J2
REAL (KIND=DBL) :: beta, tau_loc
REAL (KIND=DBL) :: T_ACC, qx_ACC, qy_ACC
REAL (KIND=DBL) :: T_VDF, qx_VDF, qy_VDF
REAL (KIND=DBL) :: l_T, l_qx, l_qy

Temp = 0.d0
Qx = 0.d0
Qy = 0.d0

T_s = 0.d0
Qx_s = 0.d0
Qy_s = 0.d0

!$OMP PARALLEL DO PRIVATE(I,M,J1,J2,beta,tau_loc, &
!$OMP T_ACC,qx_ACC,qy_ACC, T_VDF,qx_VDF,qy_VDF,l_T,l_qx,l_qy)
DO I = 1, N_TRIS
        tau_loc = TAU_R/TRIANGLES_Hmin(I)
        beta = MIN(tau_loc,TAU_THR)/tau_loc

    DO M = 1, NDOF_TRI
        T_VDF = 0.d0
        qx_VDF = 0.d0
        qy_VDF = 0.d0
        DO J2 = 1, NAZIM
        DO J1 = 1, NPOLE
            T_VDF  = T_VDF  +           VDF(M,I,J1,J2)*DOMEGA(J1,J2)
            qx_VDF = qx_VDF + CX(J1,J2)*VDF(M,I,J1,J2)*DOMEGA(J1,J2)
            qy_VDF = qy_VDF + CY(J1,J2)*VDF(M,I,J1,J2)*DOMEGA(J1,J2)
        END DO
        END DO
        T_VDF = T_VDF/Cv

        T_ACC  = UQ(           M, I)
        qx_ACC = UQ(  NDOF_TRI+M, I)
        qy_ACC = UQ(2*NDOF_TRI+M, I)

        l_T  = ( T_ACC -  T_VDF)*beta
        l_qx = (qx_ACC - qx_VDF)*beta
        l_qy = (qy_ACC - qy_VDF)*beta

        DO J2 = 1, NAZIM
        DO J1 = 1, NPOLE
            VDF(M,I,J1,J2) = VDF(M,I,J1,J2) + l_T*Cv/4.d0/PI + &
                (CX(J1,J2)*l_qx+CY(J1,J2)*l_qy)*TAU_C/TAU_N*3.d0/4.d0/PI/Vg/Vg
        END DO
        END DO

        T_s(M,I)  =  l_T +  T_VDF
        Qx_s(M,I) = l_qx + qx_VDF
        Qy_s(M,I) = l_qy + qy_VDF


        Temp(I) = Temp(I) + ( l_T +  T_VDF)*INT_NODFUNC_TRI(M,I)
        Qx(I)   =   Qx(I) + (l_qx + Qx_VDF)*INT_NODFUNC_TRI(M,I)
        Qy(I)   =   Qy(I) + (l_qy + Qy_VDF)*INT_NODFUNC_TRI(M,I)
    END DO
END DO
!$OMP END PARALLEL DO

MASS = SUM(Temp)

END SUBROUTINE

SUBROUTINE Release_PARDISO ()
IMPLICIT NONE
EXTERNAL :: PARDISO
INTEGER :: N

WRITE(*,*) 'Release PARDISO ...'
N = N_FCS*3*NDOF_FC
PHASE = -1
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,N,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)

END SUBROUTINE

END MODULE
