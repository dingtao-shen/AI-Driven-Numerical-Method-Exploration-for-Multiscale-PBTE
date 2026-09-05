MODULE VELOCITY_DISTRIBUTION_FUNCTION
USE CONSTANT
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE VELOCITY_GRID
USE BOUNDARY_CONDITION
USE INTEGRATIONS
IMPLICIT NONE

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: VDF  !perturbed VDF, DIMENSION (NDOF_TRI,N_TRIS,NPOLE,NAZIM)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: Temp, Qx, Qy !cell average value
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: T_OLD, Qx_OLD, Qy_OLD
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: T_s, Qx_s, Qy_s
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: FLUX_WALL

REAL (KIND=DBL), SAVE :: MASS

CONTAINS

SUBROUTINE Init_Velocity_Distribution_Function()
IMPLICIT NONE
INTEGER :: ios
!REAL (KIND=DBL) :: Tc, Th, x1, x2, x, Tw
CHARACTER(LEN=100) :: CHAR1

ALLOCATE(VDF(NDOF_TRI,N_TRIS,NPOLE,NAZIM))
ALLOCATE(Temp(N_TRIS),Qx(N_TRIS),Qy(N_TRIS))
ALLOCATE(T_OLD(N_TRIS),Qx_OLD(N_TRIS),Qy_OLD(N_TRIS))
ALLOCATE(T_s(NDOF_TRI,N_TRIS),Qx_s(NDOF_TRI,N_TRIS),Qy_s(NDOF_TRI,N_TRIS))
ALLOCATE(FLUX_WALL(NDOF_TRI,N_FCS_B))

WRITE(*,*)'**** Initial Velocity Distribution Function ****'
VDF = 0.d0
Temp = 0.d0
Qx = 0.d0
Qy = 0.d0

T_OLD  = 0.d0
Qx_OLD = 0.d0
Qy_OLD = 0.d0

T_s = 0.d0
Qx_s = 0.d0
Qy_s = 0.d0

FLUX_WALL = 0.d0

WRITE(CHAR1,108)'P',DEG,'T',N_TRIS,'NP',NPOLE,'NA',NAZIM
108 FORMAT(A1,I1,A1,I3,A2,I2,A2,I2)
OPEN (UNIT = 12, FILE='VDF'//TRIM(CHAR1)//'.out', STATUS='OLD', ACCESS='STREAM', FORM='UNFORMATTED',IOSTAT=ios)
IF (ios.EQ.0) THEN
WRITE(*,*)'read from file: ', 'VDF'//TRIM(CHAR1)//'.out'
READ (12) VDF
END IF
CLOSE(12)

END SUBROUTINE Init_Velocity_Distribution_Function

!************************************************************
SUBROUTINE Calculate_FLUX_WALL ()
IMPLICIT NONE
INTEGER :: I,J1,J2,M,L, BCID, TRID, LFCID
REAL (KIND=DBL) :: ss1, ss2, n1, n2, u

FLUX_WALL =0.d0

!$OMP PARALLEL DO PRIVATE (I,J1,J2,M,L,BCID,TRID,LFCID,ss1,ss2,n1,n2,u)
DO I = 1, N_FCS_B
    BCID = FACES_TAG(I,3)
    IF(BC_TYP(BCID).EQ.2)THEN
        !non-thermalizing wall
        IF (FACES_TAG(I,5).EQ.0) THEN
            TRID = FACES_TAG(I,6)
            LFCID = FACES_TAG(I,8)
        ELSE
            TRID = FACES_TAG(I,5)
            LFCID = FACES_TAG(I,7)
        END IF
        n1 = -TRIANGLES_INF(TRID,2*LFCID)
        n2 = -TRIANGLES_INF(TRID,2*LFCID+1)

        DO M = 1, NDOF_TRI
            ss1 = 0.d0
            ss2 = 0.d0
            DO J2 = 1, NAZIM
            DO J1 = 1, NPOLE
                u = CX(J1,J2)*n1+CY(J1,J2)*n2
                IF (u.LT.0.d0)THEN
                    DO L = 1, NDOF_TRI
                        ss1 = ss1 + u*DOMEGA(J1,J2)*VDF(L,TRID,J1,J2)*INT_NODFUNC_TRI_FC(TRID,LFCID,M,L)
                    END DO
                ELSE
                    ss2 = ss2 - u*DOMEGA(J1,J2)
                END IF
            END DO
            END DO
            FLUX_WALL(M,I) = -ss1/ss2
        END DO
    END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Calculate_Macro_Properties ()
IMPLICIT NONE
INTEGER :: I,L
REAL (KIND=DBL) :: ss,ss_x,ss_y

Temp = 0.d0
Qx = 0.d0
Qx = 0.d0

T_s = 0.d0
Qx_s = 0.d0
Qy_s = 0.d0

!$OMP PARALLEL DO PRIVATE (I,L,ss,ss_x,ss_y)
DO I = 1, N_TRIS
    DO L = 1, NDOF_TRI
        ss   = SUM(   VDF(L,I,:,:)*DOMEGA)
        ss_x = SUM(CX*VDF(L,I,:,:)*DOMEGA)
        ss_y = SUM(CY*VDF(L,I,:,:)*DOMEGA)

        T_s(L,I)  = ss/Cv
        Qx_s(L,I) = ss_x
        Qy_s(L,I) = ss_y

        Temp(I) = Temp(I) + ss  *INT_NODFUNC_TRI(L,I)/Cv
        Qx(I)   = Qx(I)   + ss_x*INT_NODFUNC_TRI(L,I)
        Qy(I)   = Qy(I)   + ss_y*INT_NODFUNC_TRI(L,I)
    END DO
END DO
!$OMP END PARALLEL DO

MASS = SUM(Temp)

END SUBROUTINE Calculate_Macro_Properties

SUBROUTINE Calculate_Residual_T (RESIDUAL)
IMPLICIT NONE
REAL (KIND=DBL), INTENT(OUT) :: RESIDUAL
REAL (KIND=DBL) :: ss1,ss2
INTEGER :: I

ss1 = 0.d0
ss2 = 0.d0
!$OMP PARALLEL DO PRIVATE (I) REDUCTION(+:ss1,ss2)
DO I = 1, N_TRIS
    ss1 = ss1 + (Temp(I)-T_OLD(I))*(Temp(I)-T_OLD(I))
    ss2 = ss2 + Temp(I)*Temp(I)
END DO
!$OMP END PARALLEL DO

RESIDUAL = sqrt(ss1/ss2)

T_OLD = Temp

END SUBROUTINE Calculate_Residual_T

SUBROUTINE Calculate_Residual_ALL (RESIDUAL)
IMPLICIT NONE
REAL (KIND=DBL), INTENT(OUT) :: RESIDUAL
REAL (KIND=DBL), DIMENSION(3) :: ss1,ss2
INTEGER :: I
REAL (KIND=DBL) :: q, q_old

ss1 = 0.d0
ss2 = 0.d0
!$OMP PARALLEL DO PRIVATE (I,q,q_old) REDUCTION(+:ss1,ss2)
DO I = 1, N_TRIS
    !q = SQRT(Qx(I)*Qx(I)+Qy(I)*Qy(I))
    !q_old = SQRT(Qx_OLD(I)*Qx_OLD(I)+Qy_OLD(I)*Qy_OLD(I))

    ss1(1) = ss1(1) + (Temp(I)-T_OLD(I))*(Temp(I)-T_OLD(I))
    ss1(2) = ss1(2) + (Qx(I)-Qx_old(I))*(Qx(I)-Qx_old(I))
    ss1(3) = ss1(3) + (Qy(I)-Qy_old(I))*(Qy(I)-Qy_old(I))

    ss2(1) = ss2(1) + Temp(I)*Temp(I)
    ss2(2) = ss2(2) + Qx(I)*Qx(I)
    ss2(3) = ss2(3) + Qy(I)*Qy(I)
END DO
!$OMP END PARALLEL DO

RESIDUAL = MAXVAL(SQRT(ss1/ss2))

!WRITE(*,201) SQRT(ss1/ss2)
201 FORMAT (3ES10.2)

T_OLD = Temp
Qx_OLD = Qx
Qy_OLD = Qy

END SUBROUTINE Calculate_Residual_ALL

subroutine Compress_VDF()
implicit none
CHARACTER(LEN=100) :: CHAR1

WRITE(CHAR1,108)'P',DEG,'T',N_TRIS,'NP',NPOLE,'NA',NAZIM
108 FORMAT(A1,I1,A1,I3,A2,I2,A2,I2)
OPEN (UNIT = 12, FILE='VDF'//TRIM(CHAR1)//'.out', STATUS='UNKNOWN', ACCESS='STREAM', FORM='UNFORMATTED')
WRITE(12) VDF
CLOSE(12)
WRITE(*,*)'Compress VDF'

end subroutine Compress_VDF


END MODULE VELOCITY_DISTRIBUTION_FUNCTION
