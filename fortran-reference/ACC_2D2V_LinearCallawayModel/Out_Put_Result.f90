MODULE OUT_PUT
USE CONSTANT
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE VELOCITY_GRID
USE VELOCITY_DISTRIBUTION_FUNCTION
USE ACCELERATION
IMPLICIT NONE

CONTAINS


SUBROUTINE Out_Put_Conduction()
IMPLICIT NONE
INTEGER :: Nx = 109
INTEGER :: Ny = 109
INTEGER :: I, IL, J, K,M,L,KL,JL,J1,J2
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:) :: px,py,pm
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:) :: pT, pqx, pqy, pxx, pxy, pyy
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TRIID
REAL (KIND=DBL), DIMENSION(4) :: x,y
LOGICAL :: INTRI
REAL (KIND=DBL) :: AREA, A, B, C, D, E, F, G, H, xi, eta
REAL (KIND=DBL) ::  ss,ss_x,ss_y,ss_xx,ss_xy,ss_yy

CHARACTER(LEN=100) :: FILENAME,INTG3
CHARACTER(LEN=10) :: INTG1, INTG2, INTG4

WRITE(*,*)'Output Conduction'

ALLOCATE(px(Nx),py(Ny),pm(NDOF_TRI),TRIID(Nx,Ny))
ALLOCATE(pT(Nx,Ny),pqx(Nx,Ny),pqy(Nx,Ny),pxx(Nx,Ny),pxy(Nx,Ny),pyy(Nx,Ny))

 px = 0.d0
 py = 0.d0
 pm = 0.d0
 pT = 0.d0
 pqx = 0.d0
 pqy = 0.d0
 pxx = 0.d0
 pxy = 0.d0
 pyy = 0.d0
 TRIID = 0



DO J = 1, 5
A = (J-1.d0)/8.d0
py(J) = A*A*A*(10.d0-15.d0*A+6.d0*A*A)*0.02d0
py(Ny-J+1) = 1.d0 - py(J)
END DO

DO J = 6, Ny-5
py(J) = (J-5.d0)*(1.d0-0.02d0)/100.d0 + 0.01d0
END DO
px = py

DO J = 1, Ny
    DO I = 1, Nx
        DO K = 1, N_TRIS
            x(1) = NODES(TRIANGLES_TAG(K,1),1)
            y(1) = NODES(TRIANGLES_TAG(K,1),2)
            x(2) = NODES(TRIANGLES_TAG(K,2),1)
            y(2) = NODES(TRIANGLES_TAG(K,2),2)
            x(3) = NODES(TRIANGLES_TAG(K,3),1)
            y(3) = NODES(TRIANGLES_TAG(K,3),2)
            x(4) = x(1)
            y(4) = y(1)

            INTRI = .TRUE.
            DO IL = 1,3
                AREA = 0.5*(px(I)*(y(IL)-y(IL+1))-py(J)*(x(IL)-x(IL+1))+(x(IL)*y(IL+1)-x(IL+1)*y(IL)))
                IF (AREA<-1.d-12) THEN
                    INTRI = .FALSE.
                    EXIT
                END IF
            END DO
            IF (INTRI) THEN
                TRIID(I,J) = K
                EXIT
            END IF
        END DO
    END DO
 END DO

DO J = 1, Ny
    DO I = 1, Nx
        IF (TRIID(I,J).EQ.0) THEN
            WRITE(*,*)'Warning: the point (',I,',',J,') is not in any triangle'
        ELSE
            x(1) = NODES(TRIANGLES_TAG(TRIID(I,J),1),1)
            y(1) = NODES(TRIANGLES_TAG(TRIID(I,J),1),2)
            x(2) = NODES(TRIANGLES_TAG(TRIID(I,J),2),1)
            y(2) = NODES(TRIANGLES_TAG(TRIID(I,J),2),2)
            x(3) = NODES(TRIANGLES_TAG(TRIID(I,J),3),1)
            y(3) = NODES(TRIANGLES_TAG(TRIID(I,J),3),2)
            x(4) = x(1)
            y(4) = y(1)

            A=(x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1))
            B=y(3)-y(1)
            C=x(1)-x(3)
            D=(x(3)-x(1))*y(1)-(y(3)-y(1))*x(1)

            E=(x(3)-x(1))*(y(2)-y(1))-(x(2)-x(1))*(y(3)-y(1))
            F=y(2)-y(1)
            G=x(1)-x(2)
            H=(x(2)-x(1))*y(1)-(y(2)-y(1))*x(1)

            xi = B/A*px(I)+C/A*py(J)+D/A
            eta = F/E*px(I)+G/E*py(J)+H/E

            pm = 0.d0
            DO M = 1, NDOF_TRI
                L = 1
                DO KL = 0,DEG
                    DO JL = 0, KL
                        IL = KL - JL
                        pm(M) = pm(M) + NODFUN_TRI_REF(M,L)*(xi**REAL(IL,DBL))*(eta**REAL(JL,DBL))
                        L = L + 1
                    END DO
                END DO
            END DO

            ss = 0.d0
            ss_x = 0.d0
            ss_y = 0.d0

            DO J2 = 1, NAZIM
                DO J1 = 1, NPOLE
                    DO L = 1, NDOF_TRI
                        ss   = ss   +           DOMEGA(J1,J2)*VDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ss_x = ss_x + CX(J1,J2)*DOMEGA(J1,J2)*VDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ss_y = ss_y + CY(J1,J2)*DOMEGA(J1,J2)*VDF(L,TRIID(I,J),J1,J2)*pm(L)

                    END DO
                END DO
            END DO

            pT(I,J) = ss/Cv
            pqx(I,J) = ss_x
            pqy(I,J) = ss_y

            ss_xx=0.d0
            ss_xy=0.d0
            ss_yy=0.d0
            IF (ACCFLAG.EQ.1) THEN
            DO L =1, NDOF_TRI
                ss_xx = ss_xx +(-4.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))+2.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)
                ss_xy = ss_xy +(-          UQ(4*NDOF_TRI+L,TRIID(I,J))-          UQ(5*NDOF_TRI+L,TRIID(I,J)))*pm(L)
                ss_yy = ss_yy +(+2.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))-4.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            END DO
            END IF
            pxx(I,J)=ss_xx
            pxy(I,J)=ss_xy
            pyy(I,J)=ss_yy
        END IF
    END DO
END DO

!**********************************************************
WRITE(INTG1,101)'P',DEG
101 FORMAT (A1,I1)
WRITE(INTG2,102)'T',N_TRIS
102 FORMAT (A1,I5)
WRITE(INTG3,103)'tR',TAU_R,'_tN',TAU_N
103 FORMAT (A2,ES9.2,A3,ES9.2)
IF (ACCFLAG.EQ.0) WRITE(INTG4,104) '_ CIS'
IF (ACCFLAG.EQ.1) WRITE(INTG4,104) '_GSIS'
104 FORMAT (A5)

!Field
FILENAME = '2D_'//TRIM(INTG1)//'_'//TRIM(INTG2)//'_'//TRIM(INTG3)//TRIM(INTG4)//'.dat'
OPEN(13,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE')
WRITE(13,*) 'VARIABLES="x","y","T","qx","qy","Nxx","Nxy","Nyy"'
WRITE(13,*) 'ZONE I = 109 J = 109'
DO J = 1, Ny
    DO I = 1, Nx
        WRITE(13,100) px(I),py(J),pT(I,J),pqx(I,J),pqy(I,J),pxx(I,J),pxy(I,J),pyy(I,J)
        100 FORMAT(1X,8ES16.6)
    END DO
END DO
CLOSE(13)

!analytical solution
OPEN(13,FILE='Conduction_A.dat',STATUS='UNKNOWN',ACTION='WRITE')
WRITE(13,*)'VARIABLES="x","y","T","qx","qy"'
WRITE(13,*)'ZONE I = 109 J = 109'
DO J = 1, Ny
    DO I = 1, Nx
        ss = 0.d0
        ss_x = 0.d0
        ss_y = 0.d0
        DO M = 1, 200
            ss = ss + ((-1.d0)**REAL(M+1,DBL)+1.d0)/REAL(M,DBL)*sin(REAL(M,DBL)*PI*px(I)) &
                    *sinh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
            ss_x = ss_x + ((-1.d0)**REAL(M+1,DBL)+1.d0)*cos(REAL(M,DBL)*PI*px(I)) &
                    *sinh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
            ss_y = ss_y + ((-1.d0)**REAL(M+1,DBL)+1.d0)*sin(REAL(M,DBL)*PI*px(I)) &
                    *cosh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
        END DO
        ss = ss*2.d0/PI
        ss_x = -ss_x*Cv/3.d0*TAU_R
        ss_y = -ss_y*Cv/3.d0*TAU_R

        WRITE(13,105)px(I),py(J),ss,ss_x,ss_y
        105 FORMAT(1X,5ES16.6)
    END DO
END DO
CLOSE(13)

!FILENAME = '1D_'//TRIM(INTG1)//'_'//TRIM(INTG2)//'_'//TRIM(INTG3)//TRIM(INTG4)//'.dat'
!OPEN(13,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE')
!WRITE(13,*) 'VARIABLES="x","T","qx","qy"'
!WRITE(13,*) 'ZONE I = 109'
!J=Ny/2+1
!DO I = 1, Nx
!    WRITE(13,105) px(I),pT(I,J),pqx(I,J),pqy(I,J)
!        105 FORMAT(1X,4ES16.6)
!END DO
!CLOSE(13)

DEALLOCATE(px,py,pT,pqx,pqy,pxx,pxy,pyy,TRIID)
END SUBROUTINE Out_Put_Conduction

END MODULE OUT_PUT
