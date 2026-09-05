MODULE INTEGRATIONS
!Define the matrix of integration of basis functions needed
USE CONSTANT
USE GLOBAL_VARIABLE
USE BOUNDARY_CONDITION
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE USER_DEFINE_MATH
IMPLICIT NONE

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: INT_NODFUNC_TRI_TRI ! Dimension(N_TRIS,NDOF_TRI,NDOF_TRI),
                                                                            ! integration of product of two triangle basis functions over triangle
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: INT_NODFUNC_TRI_TRI_TRI ! Dimension(N_TRIS,NDOF_TRI,NDOF_TRI,NDOF_TRI)
																				!integration of product of three triangle basis functions over triangle

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: INT_NODFUNC_TRI     ! DIMENSION(N_TRIS,NDOF_TRI)
                                                                            ! integration of triangle basis over triangle
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: INT_NODFUNC_TRI_TRI_X ! DIMENSION(N_TRIS,NDOF_TRI,NDOF_TRI)
                                                                              ! integration of product of one basis function and partial x derivative of one basis function over triangle
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: INT_NODFUNC_TRI_TRI_Y ! DIMENSION(N_TRIS,NDOF_TRI,NDOF_TRI)
                                                                            ! integration of product of one basis function and partial y derivative of one basis function over triangle
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: INT_NODFUNC_TRI_FC  !Dimension(N_TRIS,3,NDOF_TRI, NDOF_TRI)
                                                                              !integration of product of two triangle basis functions over each edge
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: INT_NODFUNC_TRI_FC_FC ! DIMENSION(N_TRIS,3,NDOF_TRI, NDOF_FC)
                                                                             !integration of product of one triangle basis function and one face basis function over each edge
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: INT_NODFUNC_FC_FC    !Dimension (N_FCS,NDOF_FC,NDOF_FC)
                                                                             !Integration of product of two face basis functions over face
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: INT_NODFUNC_FC         !DIMENSION(N_FCS,NDOF_FC)
                                                                             !Integration of face basis over face
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: INT_NODFUNC_TRI_TRI_FC  !DIMESION(N_TRIS,3,NDOF_TRI,NDOF_TRI)
                                                                                  !inregration of production of interior and exterior triangle basis functions over each edge

!***********************************************************
!numerical quadrature
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: QUA_ABS_X, QUA_ABS_Y, QUA_WEI !triangle quadrature in reference triangle
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: NODFUN_QUA_P ! nodal function evaluated at the quadrature point in reference triangle

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: QUA_ABS_LINE, QUA_WEI_LINE !quadrature in reference segment [-1,1]
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: NODFUN_QUA_P_LINE !nodal function evaluated at the line quadrature point
!***********************************************************

CONTAINS

	SUBROUTINE Init_Integration ()
	IMPLICIT NONE
!	REAL (KIND=DBL) :: area
	INTEGER :: I, M, L, N, KL,KM,KN,IL,IM,IN,JL,JM,JN,J1,J2, FCID,J,TRIID, BCID, PAIR, BCID_P
	REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: INT_ELMENT, INT_ELMENT_X, INT_ELMENT_Y
	REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:) :: INT_ELMENT_E, INT_ELMENT_P
	REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:) :: ABSC,WEIG,pm,pl
	REAL (KIND=DBL) :: ss, ss1,ss2,ss3,x1,x2,x3,y1,y2,y3,len,xp,yp,A,B,C,D,E,F,G,H,xi,eta
	REAL (KIND=DBL), DIMENSION(4) :: x,y
	REAL*8 :: FAC_DIV

	ALLOCATE (INT_NODFUNC_TRI(NDOF_TRI,N_TRIS))
	ALLOCATE (INT_NODFUNC_TRI_TRI(NDOF_TRI,NDOF_TRI,N_TRIS))
	ALLOCATE (INT_NODFUNC_TRI_TRI_X(N_TRIS,NDOF_TRI,NDOF_TRI))
	ALLOCATE (INT_NODFUNC_TRI_TRI_Y(N_TRIS,NDOF_TRI,NDOF_TRI))
	ALLOCATE (INT_NODFUNC_TRI_TRI_TRI(NDOF_TRI,NDOF_TRI,NDOF_TRI,N_TRIS))

	ALLOCATE (INT_NODFUNC_TRI_FC(N_TRIS,3,NDOF_TRI,NDOF_TRI))
	ALLOCATE (INT_NODFUNC_TRI_FC_FC(NDOF_TRI,NDOF_FC,N_TRIS,3))
	ALLOCATE (INT_NODFUNC_FC_FC(NDOF_FC,NDOF_FC,N_FCS))
	ALLOCATE (INT_NODFUNC_FC(NDOF_FC,N_FCS))

	ALLOCATE (INT_ELMENT(NDOF_TRI,NDOF_TRI),INT_ELMENT_E(3,NDOF_TRI,NDOF_TRI))
	ALLOCATE (INT_ELMENT_X(NDOF_TRI,NDOF_TRI),INT_ELMENT_Y(NDOF_TRI,NDOF_TRI))
	ALLOCATE (INT_ELMENT_P(NDOF_TRI,NDOF_TRI,NDOF_TRI))

	!calculate INT_NODFUNC_TRI
	INT_ELMENT = 0.d0
	INT_NODFUNC_TRI = 0.d0
	M = 1
	DO KM = 0,DEG
		DO JM = 0, KM
			IM = KM-JM
			!INT_ELMENT(M,1) = REAL(FAC(IM),DBL)*REAL(FAC(JM),DBL)/REAL(FAC(IM+JM+2))
			INT_ELMENT(M,1) = FAC_DIV(1,JM,IM+1,IM+JM+2)
			M = M + 1
		END DO
	END DO

	DO M = 1, NDOF_TRI
		ss = 0.d0
		DO IM = 1, NDOF_TRI
			ss = ss + NODFUN_TRI_REF(M,IM)*INT_ELMENT(IM,1)
		END DO
		INT_NODFUNC_TRI(M,:) = ss
	END DO

	DO I = 1, N_TRIS
		INT_NODFUNC_TRI(:,I) = INT_NODFUNC_TRI(:,I)*2.d0*TRIANGLES_INF(I,1)
	END DO

	!calculate INT_NODFUNC_TRI_TRI
	INT_ELMENT = 0.d0
	INT_NODFUNC_TRI_TRI = 0.d0
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM-JM

			L = 1
			DO KL = 0, DEG
				DO JL = 0, KL
					IL = KL - JL
					!INT_ELMENT(M,L) = REAL(FAC(IL+IM),DBL)*REAL(FAC(JL+JM),DBL)/REAL(FAC(IL+IM+JL+JM+2),DBL)
					INT_ELMENT(M,L) = FAC_DIV(1,JL+JM,IL+IM+1,IL+IM+JL+JM+2)
					L = L + 1
				END DO
			END DO
			M = M + 1
		END DO
	END DO

	DO M = 1, NDOF_TRI
		DO L = 1, NDOF_TRI
			ss=0.0
			DO IL = 1, NDOF_TRI
				DO IM = 1, NDOF_TRI
					ss = ss + NODFUN_TRI_REF (L,IL)*NODFUN_TRI_REF(M,IM)*INT_ELMENT(IM,IL)
				END DO
			END DO
			INT_NODFUNC_TRI_TRI(M,L,:) = ss
		END DO
	END DO

	DO I =1, N_TRIS
		INT_NODFUNC_TRI_TRI(:,:,I) = INT_NODFUNC_TRI_TRI(:,:,I)*2.0*TRIANGLES_INF(I,1)
	END DO

	!calculate INT_NODFUNC_TRI_TRI_TRI
	INT_ELMENT_P = 0.d0
	INT_NODFUNC_TRI_TRI_TRI = 0.d0
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM - JM

			L = 1
			DO KL = 0, DEG
				DO JL = 0, KL
					IL = KL - JL

					N = 1
					DO KN = 0, DEG
						DO JN = 0, KN
							IN = KN - JN
							!INT_ELMENT_P(M,L,N) = REAL(FAC(IL+IM+IN),DBL)*REAL(FAC(JL+JM+JN),DBL) &
							!             /REAL(FAC(IL+IM+IN+JL+JM+JN+2),DBL)
							INT_ELMENT_P(M,L,N) = FAC_DIV(1,JL+JM+JN,IL+IM+IN+1,IL+IM+IN+JL+JM+JN+2)
							N = N + 1
						END DO
					END DO

					L = L + 1
				END DO
			END DO

			M = M + 1
		END DO
	END DO


	DO M = 1, NDOF_TRI
		DO L = 1, NDOF_TRI
			DO N = 1, NDOF_TRI
				ss = 0.d0
				DO IN = 1, NDOF_TRI
					DO IL = 1, NDOF_TRI
						DO IM = 1, NDOF_TRI
							ss = ss + NODFUN_TRI_REF(N,IN)*NODFUN_TRI_REF(L,IL) &
							          *NODFUN_TRI_REF(M,IM)*INT_ELMENT_P(IM,IL,IN)
						END DO
					END DO
				END DO
				INT_NODFUNC_TRI_TRI_TRI(M,L,N,:) = ss
			END DO
		END DO
	END DO

	DO I = 1, N_TRIS
		INT_NODFUNC_TRI_TRI_TRI(:,:,:,I) = INT_NODFUNC_TRI_TRI_TRI(:,:,:,I)*2.d0*TRIANGLES_INF(I,1)
	END DO


	!calculate INT_NODFUNC_TRI_TRI_X
	!calculate INT_NODFUNC_TRI_TRI_Y
	INT_ELMENT = 0.d0
	INT_ELMENT_X = 0.d0
	INT_NODFUNC_TRI_TRI_X = 0.d0
	INT_ELMENT_Y = 0.d0
	INT_NODFUNC_TRI_TRI_Y = 0.d0
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM-JM

			IF (IM>0) THEN
				L = 1
				DO KL = 0, DEG
					DO JL = 0, KL
						IL = KL - JL
						!INT_ELMENT(M,L) = REAL(IM,DBL)*REAL(FAC(IL+IM-1),DBL)*REAL(FAC(JL+JM),DBL)/REAL(FAC(IL+IM-1+JL+JM+2))
						INT_ELMENT(M,L) = REAL(IM,DBL)*FAC_DIV(1,JL+JM,IL+IM,IL+IM-1+JL+JM+2)
						L = L + 1
					END DO
				END DO
			END IF
			M = M + 1
		END DO
	END DO

	DO M = 1, NDOF_TRI
		DO L = 1, NDOF_TRI
			ss=0.0
			DO IL = 1, NDOF_TRI
				DO IM = 1, NDOF_TRI
					ss = ss + NODFUN_TRI_REF (L,IL)*NODFUN_TRI_REF(M,IM)*INT_ELMENT(IM,IL)
				END DO
			END DO
			INT_ELMENT_X(M,L) = ss
		END DO
	END DO

	INT_ELMENT = 0.d0
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM-JM

			IF (JM>0) THEN
				L = 1
				DO KL = 0, DEG
					DO JL = 0, KL
						IL = KL - JL
						!INT_ELMENT(M,L) = REAL(JM,DBL)*REAL(FAC(IL+IM),DBL)*REAL(FAC(JL+JM-1),DBL)/REAL(FAC(IL+IM-1+JL+JM+2))
						INT_ELMENT(M,L) = REAL(JM,DBL)*FAC_DIV(1,JL+JM-1,IL+IM+1,IL+IM-1+JL+JM+2)
						L = L + 1
					END DO
				END DO
			END IF
			M = M + 1
		END DO
	END DO

	DO M = 1, NDOF_TRI
		DO L = 1, NDOF_TRI
			ss=0.0
			DO IL = 1, NDOF_TRI
				DO IM = 1, NDOF_TRI
					ss = ss + NODFUN_TRI_REF (L,IL)*NODFUN_TRI_REF(M,IM)*INT_ELMENT(IM,IL)
				END DO
			END DO
			INT_ELMENT_Y(M,L) = ss
		END DO
	END DO

	DO I =1, N_TRIS
		x1 = NODES (TRIANGLES_TAG(I,1),1)
		x2 = NODES (TRIANGLES_TAG(I,2),1)
		x3 = NODES (TRIANGLES_TAG(I,3),1)
		y1 = NODES (TRIANGLES_TAG(I,1),2)
		y2 = NODES (TRIANGLES_TAG(I,2),2)
		y3 = NODES (TRIANGLES_TAG(I,3),2)

		ss = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

		INT_NODFUNC_TRI_TRI_X (I,:,:) = ((y3-y1)/ss*INT_ELMENT_X+(y1-y2)/ss*INT_ELMENT_Y)*2.0*TRIANGLES_INF(I,1)
		INT_NODFUNC_TRI_TRI_Y (I,:,:) = ((x1-x3)/ss*INT_ELMENT_X+(x2-x1)/ss*INT_ELMENT_Y)*2.0*TRIANGLES_INF(I,1)
	END DO

	!calculate INT_NODFUNC_TRI_FC
	INT_NODFUNC_TRI_FC = 0.d0
	INT_ELMENT_E = 0.d0
	!first line: y = 0 in reference triangle
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM - JM
			IF (JM.EQ.0)THEN
				L = 1
				DO KL = 0, DEG
					DO JL = 0, KL
						IL = KL - JL
						IF(JL.EQ.0)THEN
							INT_ELMENT_E(1,M,L) = 1.d0/REAL(IM+IL+1,DBL)
						END IF
						L = L + 1
					END DO
				END DO
			END IF
			M = M + 1
		END DO
	END DO
	!second line: y=-x in reference triangle
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM-JM

			L = 1
			DO KL = 0, DEG
				DO JL = 0, KL
					IL = KL - JL
					!INT_ELMENT_E(2,M,L) = sqrt(2.d0)*REAL(FAC(IM+IL),DBL)*REAL(FAC(JM+JL),DBL)/&
					!     REAL(FAC(IM+IL+JM+JL+1),DBL)
					INT_ELMENT_E(2,M,L) = sqrt(2.d0)*FAC_DIV(1,JM+JL,IM+IL+1,IM+IL+JM+JL+1)
					L = L + 1
				END DO
			END DO
			M = M + 1
		END DO
	END DO
	!third line : x=0 in reference triangle
	M = 1
	DO KM = 0, DEG
		DO JM = 0, KM
			IM = KM-JM
			IF (IM.EQ.0) THEN

				L = 1
				DO KL = 0, DEG
					DO JL = 0, KL
						IL = KL -JL
						IF (IL.EQ.0)THEN
							INT_ELMENT_E(3,M,L) = 1.d0/REAL(JM+JL+1,DBL)
						END IF
						L= L + 1
					END DO
				END DO
			END IF
			M = M + 1
		END DO
	END DO

	DO I = 1, N_TRIS
		x1 = NODES (TRIANGLES_TAG(I,1),1)
		x2 = NODES (TRIANGLES_TAG(I,2),1)
		x3 = NODES (TRIANGLES_TAG(I,3),1)
		y1 = NODES (TRIANGLES_TAG(I,1),2)
		y2 = NODES (TRIANGLES_TAG(I,2),2)
		y3 = NODES (TRIANGLES_TAG(I,3),2)

		DO M = 1, NDOF_TRI
			DO L = 1, NDOF_TRI
				ss1 = 0.d0
				ss2 = 0.d0
				ss3 = 0.d0
				DO IM = 1, NDOF_TRI
					DO IL = 1, NDOF_TRI
						ss1 = ss1 + NODFUN_TRI_REF(M,IM)*NODFUN_TRI_REF(L,IL)*INT_ELMENT_E(1,IM,IL)
						ss2 = ss2 + NODFUN_TRI_REF(M,IM)*NODFUN_TRI_REF(L,IL)*INT_ELMENT_E(2,IM,IL)
						ss3 = ss3 + NODFUN_TRI_REF(M,IM)*NODFUN_TRI_REF(L,IL)*INT_ELMENT_E(3,IM,IL)
					END DO
				END DO
				INT_NODFUNC_TRI_FC(I,1,M,L) = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*ss1
				INT_NODFUNC_TRI_FC(I,2,M,L) = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))/sqrt(2.d0)*ss2
				INT_NODFUNC_TRI_FC(I,3,M,L) = sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))*ss3
			END DO
		END DO
	END DO

!	!calculate INT_NODFUNC_TRI_FC_FC
!	INT_ELMENT_E = 0.d0
!	INT_NODFUNC_TRI_FC_FC = 0.d0
!	!first line in reference triangle: y=0
!	M = 1
!	DO KM = 0, DEG
!		DO JM = 0, KM
!			IM = KM-JM
!			IF (JM.EQ.0)THEN
!				DO L = 0, DEG
!					SELECT CASE (IM)
!					CASE (0)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(1,M,L+1) = 2.d0/REAL(L+1,DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(1,M,L+1) = 0.d0
!					CASE (1)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(1,M,L+1) = 1.d0/REAL(L+1,DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(1,M,L+1) = 1.d0/REAL(L+2,DBL)
!					CASE (2)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(1,M,L+1) = REAL(4+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(1,M,L+1) = REAL(3+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!					CASE (3)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(1,M,L+1) = REAL(24+L*(34+L*(15+2*L)),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(1,M,L+1) = REAL(21+L*(34+L*(15+2*L)),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!					CASE (4)
!						IF(MOD(L,2).EQ.0) INT_ELMENT_E(1,M,L+1) = REAL(48+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!						IF(MOD(L,2).EQ.1) INT_ELMENT_E(1,M,L+1) = REAL(45+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!					END SELECT
!				END DO
!			END IF
!
!			M = M + 1
!		END DO
!	END DO

	!second line: y=1-x in reference triangle
!	DO L = 0, DEG
!		IF (MOD(L,2).EQ.0)THEN
!			INT_ELMENT_E (2,1,L+1) = 2.d0/REAL(L+1,DBL)
!			INT_ELMENT_E (2,2,L+1) = 1.d0/REAL(L+1,DBL)
!			INT_ELMENT_E (2,3,L+1) = 1.d0/REAL(L+1,DBL)
!			IF (DEG.GT.1)THEN
!				INT_ELMENT_E (2,4,L+1) = REAL(4+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!				INT_ELMENT_E (2,5,L+1) = 1.d0/REAL((L+1)*(L+3),DBL)
!				INT_ELMENT_E (2,6,L+1) = REAL(4+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!			END IF
!			IF (DEG.GT.2)THEN
!				INT_ELMENT_E (2,7,L+1) = REAL(3+(5+2*L)*(9+2*L*(5+L)),DBL)/REAL(4*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,8,L+1) = REAL(8+L+L*(5+L),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,9,L+1) = REAL(8+L*(6+L),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,10,L+1) = REAL(24+L*(34+L*(15+2*L)),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!			END IF
!			IF (DEG.GT.3) THEN
!				INT_ELMENT_E (2,11,L+1) = REAL(48+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,12,L+1) = REAL(8+L*(6+L),DBL)/REAL(2*(L+1)*(L+2)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,13,L+1) = 1.d0/REAL((L+1)*(L+3)*(L+5),DBL)
!				INT_ELMENT_E (2,14,L+1) = REAL(8+L*(6+L),DBL)/REAL(2*(L+1)*(L+2)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,15,L+1) = REAL(48+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!			END IF
!		END IF
!		IF (MOD(L,2).EQ.1)THEN
!			INT_ELMENT_E (2,1,L+1) = 0.d0
!			INT_ELMENT_E (2,2,L+1) = -1.d0/REAL(L+2,DBL)
!			INT_ELMENT_E (2,3,L+1) = 1.d0/REAL(L+2,DBL)
!			IF (DEG.GT.1)THEN
!				INT_ELMENT_E (2,4,L+1) = -REAL(3+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!				INT_ELMENT_E (2,5,L+1) = 0.d0
!				INT_ELMENT_E (2,6,L+1) = REAL(3+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!			END IF
!			IF (DEG.GT.2)THEN
!				INT_ELMENT_E (2,7,L+1) = REAL(3-(5+2*L)*(9+2*L*(5+L)),DBL)/REAL(4*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,8,L+1) = REAL(L-3-L*(5+L),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,9,L+1) = REAL(3+L*(4+L),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!				INT_ELMENT_E (2,10,L+1) = REAL(21+L*(34+L*(15+2*L)),DBL)/REAL(2*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!			END IF
!			IF (DEG.GT.3) THEN
!				INT_ELMENT_E (2,11,L+1) = REAL(-(45+L*(6+L)*(14+L*(6+L))),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,12,L+1) = -REAL(5+L*(6+L),DBL)/REAL(2*(L+1)*(L+2)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,13,L+1) = 0.d0
!				INT_ELMENT_E (2,14,L+1) = REAL(5+L*(6+L),DBL)/REAL(2*(L+1)*(L+2)*(L+4)*(L+5),DBL)
!				INT_ELMENT_E (2,15,L+1) = REAL(45+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!			END IF
!		END IF
!	END DO

	!third line : x=0 in reference triangle
!	M = 1
!	DO KM = 0,DEG
!		DO JM = 0,KM
!			IM = KM-JM
!			IF (IM.EQ.0)THEN
!				DO L = 0, DEG
!					SELECT CASE (JM)
!					CASE (0)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(3,M,L+1) = 2.d0/REAL(L+1,DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(3,M,L+1) = 0.d0
!					CASE (1)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(3,M,L+1) = 1.d0/REAL(L+1,DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(3,M,L+1) = -1.d0/REAL(L+2,DBL)
!					CASE (2)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(3,M,L+1) = REAL(4+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(3,M,L+1) = -REAL(3+L*(4+L),DBL)/REAL((L+1)*(L+2)*(L+3),DBL)
!					CASE (3)
!						IF (MOD(L,2).EQ.0) INT_ELMENT_E(3,M,L+1) = REAL(3+(5+2*L)*(9+2*L*(5+L)),DBL)/REAL(4*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!						IF (MOD(L,2).EQ.1) INT_ELMENT_E(3,M,L+1) = REAL(3-(5+2*L)*(9+2*L*(5+L)),DBL)/REAL(4*(L+1)*(L+2)*(L+3)*(L+4),DBL)
!					CASE (4)
!						IF(MOD(L,2).EQ.0) INT_ELMENT_E(3,M,L+1) = REAL(48+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!						IF(MOD(L,2).EQ.1) INT_ELMENT_E(3,M,L+1) = -REAL(45+L*(6+L)*(14+L*(6+L)),DBL)/REAL((L+1)*(L+2)*(L+3)*(L+4)*(L+5),DBL)
!					END SELECT
!				END DO
!			END IF
!			M = M + 1
!		END DO
!	END DO
!
!	DO I = 1, N_TRIS
!		x1 = NODES (TRIANGLES_TAG(I,1),1)
!		x2 = NODES (TRIANGLES_TAG(I,2),1)
!		x3 = NODES (TRIANGLES_TAG(I,3),1)
!		y1 = NODES (TRIANGLES_TAG(I,1),2)
!		y2 = NODES (TRIANGLES_TAG(I,2),2)
!		y3 = NODES (TRIANGLES_TAG(I,3),2)
!
!		DO M = 1, NDOF_TRI
!			DO L = 1, NDOF_FC
!				ss1 = 0.d0
!				ss2 = 0.d0
!				ss3 = 0.d0
!				DO IM = 1, NDOF_TRI
!					DO IL = 1, NDOF_FC
!						ss1 = ss1 + NODFUN_TRI_REF(M,IM)*NODFUN_FC_REF(L,IL)*INT_ELMENT_E(1,IM,IL)
!						ss2 = ss2 + NODFUN_TRI_REF(M,IM)*NODFUN_FC_REF(L,IL)*INT_ELMENT_E(2,IM,IL)
!						ss3 = ss3 + NODFUN_TRI_REF(M,IM)*NODFUN_FC_REF(L,IL)*INT_ELMENT_E(3,IM,IL)
!					END DO
!				END DO
!				INT_NODFUNC_TRI_FC_FC(M,L,I,1) = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*0.5*ss1
!				INT_NODFUNC_TRI_FC_FC(M,L,I,2) = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*0.5*ss2
!				INT_NODFUNC_TRI_FC_FC(M,L,I,3) = sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))*0.5*ss3
!			END DO
!
!			DO IL = 1,3
!				FCID = TRIANGLES_TAG(I,3+IL)
!				IF (FACES_TAG(FCID,1).NE.TRIANGLES_TAG(I,IL)) THEN
!					DO L = 1, NDOF_FC/2
!						ss = INT_NODFUNC_TRI_FC_FC(M,L,I,IL)
!						INT_NODFUNC_TRI_FC_FC(M,L,I,IL) = INT_NODFUNC_TRI_FC_FC(M,NDOF_FC+1-L,I,IL)
!						INT_NODFUNC_TRI_FC_FC(M,NDOF_FC+1-L,I,IL) = ss
!					END DO
!				END IF
!			END DO
!		END DO
!	END DO
!************************************************************************************************
	!numerical integration
	INT_NODFUNC_TRI_FC_FC = 0.d0
	ALLOCATE (ABSC(NP_FC),WEIG(NP_FC),pm(NDOF_TRI),pl(NDOF_FC))
	ABSC = 0.d0
	WEIG = 0.d0
	pm = 0.d0
	pl = 0.d0

	DO I = 1, N_TRIS
		!triangle
		x(1) = NODES(TRIANGLES_TAG(I,1),1)
		y(1) = NODES(TRIANGLES_TAG(I,1),2)
		x(2) = NODES(TRIANGLES_TAG(I,2),1)
		y(2) = NODES(TRIANGLES_TAG(I,2),2)
		x(3) = NODES(TRIANGLES_TAG(I,3),1)
		y(3) = NODES(TRIANGLES_TAG(I,3),2)
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

		DO J2 = 1, 3
			FCID = TRIANGLES_TAG(I,3+J2)

			len = FACES_INF(FCID)
			x1 = NODES (FACES_TAG(FCID,1),1)
			y1 = NODES (FACES_TAG(FCID,1),2)
			x2 = NODES (FACES_TAG(FCID,2),1)
			y2 = NODES (FACES_TAG(FCID,2),2)

			CALL GaussLegendre (NP_FC,0.d0,len,ABSC,WEIG)

			DO J1 = 1, NP_FC
				!triangle in reference
				xp = (x2-x1)/len*ABSC(J1) + x1
				yp = (y2-y1)/len*ABSC(J1) + y1
				xi = B/A*xp+C/A*yp+D/A
				eta = F/E*xp+G/E*yp+H/E

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

				!face in reference
				xi = 2.d0*ABSC(J1)/len - 1.d0

				pl = 0.d0
				DO L = 1, NDOF_FC
					DO IL = 0, DEG
						pl(L) = pl(L) + NODFUN_FC_REF(L,IL+1)*(xi**REAL(IL,DBL))
					END DO
				END DO

				DO M = 1, NDOF_TRI
					DO L = 1, NDOF_FC
						INT_NODFUNC_TRI_FC_FC(M,L,I,J2) = INT_NODFUNC_TRI_FC_FC(M,L,I,J2) + pm(M)*pl(L)*WEIG(J1)
					END DO
				END DO
			END DO
		END DO
	END DO

	DEALLOCATE(ABSC,WEIG,pm,pl)
	!**************************************************************************************************

	!calculate INT_NODFUNC_FC_FC
	INT_NODFUNC_FC_FC = 0.d0
	INT_ELMENT = 0.d0
	DO M = 0, DEG
		DO L = 0, DEG
			IF (MOD(M+L,2).EQ.0) INT_ELMENT(M+1,L+1) = 2.d0/REAL(M+L+1,DBL)
		END DO
	END DO

	DO I = 1, N_FCS
		DO M = 1, NDOF_FC
			DO L = 1, NDOF_FC
				ss = 0.d0
				DO IM = 1, NDOF_FC
					DO IL = 1, NDOF_FC
						ss = ss + NODFUN_FC_REF(M,IM)*NODFUN_FC_REF(L,IL)*INT_ELMENT(IM,IL)
					END DO
				END DO
				INT_NODFUNC_FC_FC(M,L,I) = FACES_INF(I)*0.5*ss
			END DO
		END DO
	END DO

	!calculate INT_NODFUNC_FC
	DO I = 1, N_FCS
		DO M = 1, NDOF_FC
			ss = 0.d0
			DO IM = 1, NDOF_FC
				ss = ss + NODFUN_FC_REF(M,IM)*INT_ELMENT(IM,1)
			END DO
			INT_NODFUNC_FC(M,I) = FACES_INF(I)*0.5*ss
		END DO
	END DO

	!*****************************************************************************************
	!Numerical quadrature
	ALLOCATE(QUA_ABS_X(NP_TRI),QUA_ABS_Y(NP_TRI),QUA_WEI(NP_TRI))
	ALLOCATE(NODFUN_QUA_P(NDOF_TRI,NP_TRI))

	QUA_WEI = 0.d0
	QUA_ABS_X = 0.d0
	QUA_ABS_Y = 0.d0
	CALL TRI_QUADRATURE (QUA_ABS_X,QUA_ABS_Y,QUA_WEI,NP_TRI)

	NODFUN_QUA_P = 0
	DO M = 1, NDOF_TRI
		DO J = 1, NP_TRI
			L = 1
			DO KL = 0, DEG
				DO JL = 0, KL
					IL = KL - JL
					NODFUN_QUA_P(M,J) = NODFUN_QUA_P(M,J) + &
					           NODFUN_TRI_REF(M,L)*(QUA_ABS_X(J)**REAL(IL,DBL))*(QUA_ABS_Y(J)**REAL(JL,DBL))
					L = L + 1
				END DO
			END DO
		END DO
	END DO

	!****************************************************************************************

    DEALLOCATE (INT_ELMENT,INT_ELMENT_X,INT_ELMENT_Y)
	DEALLOCATE (INT_ELMENT_E,INT_ELMENT_P)

	!*************************************************************
	!calculate INT_NODFUNC_TRI_TRI_FC
	ALLOCATE (INT_NODFUNC_TRI_TRI_FC(N_TRIS,3,NDOF_TRI,NDOF_TRI))
	ALLOCATE (ABSC(NP_FC),WEIG(NP_FC),pm(NDOF_TRI),pl(NDOF_TRI))
	INT_NODFUNC_TRI_TRI_FC = 0.d0
	ABSC = 0.d0
	WEIG = 0.d0

	DO I = 1, N_TRIS
		DO J2 = 1, 3
			FCID = TRIANGLES_TAG(I,3+J2)
			BCID = FACES_TAG(FCID,3)
			
			len = FACES_INF(FCID)
			x1 = NODES (FACES_TAG(FCID,1),1)
			y1 = NODES (FACES_TAG(FCID,1),2)
			x2 = NODES (FACES_TAG(FCID,2),1)
			y2 = NODES (FACES_TAG(FCID,2),2)

			CALL GaussLegendre (NP_FC,0.d0,len,ABSC,WEIG)

			IF (BCID.EQ.0) THEN
				
				!interior
				IF (FACES_TAG(FCID,5).EQ.I) TRIID = FACES_TAG(FCID,6) !neighbor triangle
				IF (FACES_TAG(FCID,6).EQ.I) TRIID = FACES_TAG(FCID,5)

				DO J1 = 1, NP_FC
					xp = (x2-x1)/len*ABSC(J1) + x1
					yp = (y2-y1)/len*ABSC(J1) + y1
					!****************************************************
					!current triangle
					x(1) = NODES(TRIANGLES_TAG(I,1),1)
					y(1) = NODES(TRIANGLES_TAG(I,1),2)
					x(2) = NODES(TRIANGLES_TAG(I,2),1)
					y(2) = NODES(TRIANGLES_TAG(I,2),2)
					x(3) = NODES(TRIANGLES_TAG(I,3),1)
					y(3) = NODES(TRIANGLES_TAG(I,3),2)
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

					xi = B/A*xp+C/A*yp+D/A
					eta = F/E*xp+G/E*yp+H/E

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
					!****************************************************
					!neighbor triangle
					x(1) = NODES(TRIANGLES_TAG(TRIID,1),1)
					y(1) = NODES(TRIANGLES_TAG(TRIID,1),2)
					x(2) = NODES(TRIANGLES_TAG(TRIID,2),1)
					y(2) = NODES(TRIANGLES_TAG(TRIID,2),2)
					x(3) = NODES(TRIANGLES_TAG(TRIID,3),1)
					y(3) = NODES(TRIANGLES_TAG(TRIID,3),2)
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

					xi = B/A*xp+C/A*yp+D/A
					eta = F/E*xp+G/E*yp+H/E

					pl = 0.d0
					DO M = 1, NDOF_TRI
						L = 1
						DO KL = 0,DEG
							DO JL = 0, KL
								IL = KL - JL
								pl(M) = pl(M) + NODFUN_TRI_REF(M,L)*(xi**REAL(IL,DBL))*(eta**REAL(JL,DBL))
								L = L + 1
							END DO
						END DO
					END DO
					!**************************
					DO M = 1, NDOF_TRI
						DO L = 1, NDOF_TRI
							INT_NODFUNC_TRI_TRI_FC(I,J2,M,L) = INT_NODFUNC_TRI_TRI_FC(I,J2,M,L) + &
								pm(M)*pl(L)*WEIG(J1)
						END DO
					END DO
				END DO
			ELSE
				IF (BC_TYP(BCID).EQ.3) THEN
					!periodic BC
					PAIR = FACES_TAG(FCID,4)
					BCID_P = FACES_TAG(PAIR,3)
					!interior
					IF (FACES_TAG(PAIR,5).EQ.0) TRIID = FACES_TAG(PAIR,6) !neighbor triangle
					IF (FACES_TAG(PAIR,6).EQ.0) TRIID = FACES_TAG(PAIR,5)

					DO J1 = 1, NP_FC
						xp = (x2-x1)/len*ABSC(J1) + x1
						yp = (y2-y1)/len*ABSC(J1) + y1
						!****************************************************
						!current triangle
						x(1) = NODES(TRIANGLES_TAG(I,1),1)
						y(1) = NODES(TRIANGLES_TAG(I,1),2)
						x(2) = NODES(TRIANGLES_TAG(I,2),1)
						y(2) = NODES(TRIANGLES_TAG(I,2),2)
						x(3) = NODES(TRIANGLES_TAG(I,3),1)
						y(3) = NODES(TRIANGLES_TAG(I,3),2)
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

						xi = B/A*xp+C/A*yp+D/A
						eta = F/E*xp+G/E*yp+H/E

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
					!****************************************************
						!neighbor triangle
						x(1) = NODES(TRIANGLES_TAG(TRIID,1),1) + BC_XOFF(BCID_P)
						y(1) = NODES(TRIANGLES_TAG(TRIID,1),2) + BC_YOFF(BCID_P)
						x(2) = NODES(TRIANGLES_TAG(TRIID,2),1) + BC_XOFF(BCID_P)
						y(2) = NODES(TRIANGLES_TAG(TRIID,2),2) + BC_YOFF(BCID_P)
						x(3) = NODES(TRIANGLES_TAG(TRIID,3),1) + BC_XOFF(BCID_P)
						y(3) = NODES(TRIANGLES_TAG(TRIID,3),2) + BC_YOFF(BCID_P)
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

						xi = B/A*xp+C/A*yp+D/A
						eta = F/E*xp+G/E*yp+H/E

						pl = 0.d0
						DO M = 1, NDOF_TRI
							L = 1
							DO KL = 0,DEG
								DO JL = 0, KL
									IL = KL - JL
									pl(M) = pl(M) + NODFUN_TRI_REF(M,L)*(xi**REAL(IL,DBL))*(eta**REAL(JL,DBL))
									L = L + 1
								END DO
							END DO
						END DO
						!**************************
						DO M = 1, NDOF_TRI
							DO L = 1, NDOF_TRI
								INT_NODFUNC_TRI_TRI_FC(I,J2,M,L) = INT_NODFUNC_TRI_TRI_FC(I,J2,M,L) + &
									pm(M)*pl(L)*WEIG(J1)
							END DO
						END DO
					END DO
				ELSE
					!boundary exclude period
					DO J1 = 1, NP_FC
						xp = (x2-x1)/len*ABSC(J1) + x1
						yp = (y2-y1)/len*ABSC(J1) + y1
						!****************************************************
						!current triangle
						x(1) = NODES(TRIANGLES_TAG(I,1),1)
						y(1) = NODES(TRIANGLES_TAG(I,1),2)
						x(2) = NODES(TRIANGLES_TAG(I,2),1)
						y(2) = NODES(TRIANGLES_TAG(I,2),2)
						x(3) = NODES(TRIANGLES_TAG(I,3),1)
						y(3) = NODES(TRIANGLES_TAG(I,3),2)
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

						xi = B/A*xp+C/A*yp+D/A
						eta = F/E*xp+G/E*yp+H/E

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
						!****************************************************
						DO M = 1, NDOF_TRI
							INT_NODFUNC_TRI_TRI_FC(I,J2,M,:) = INT_NODFUNC_TRI_TRI_FC(I,J2,M,:) + &
								pm(M)*WEIG(J1)
						END DO
					END DO
				END IF
			END IF
		END DO
	END DO

	!***************************************************************
	!IF (TCASE.EQ.5) THEN
	!NP_FC = 2*DEG + 1
	!ALLOCATE (QUA_ABS_LINE(NP_FC),QUA_WEI_LINE(NP_FC))
	!CALL GaussLegendre(NP_FC,-1.d0,1.d0,QUA_ABS_LINE,QUA_WEI_LINE)

	!ALLOCATE (NODFUN_QUA_P_LINE(NP_FC,NDOF_TRI,N_FCS_B))

	!DO I = 1, N_FCS_B
	!IF ((FACES_TAG(I,3).EQ.1).OR.(FACES_TAG(I,3).EQ.2)) THEN
	!	IF (FACES_TAG(I,5).EQ.0) TRIID = FACES_TAG(I,6)
	!	IF (FACES_TAG(I,6).EQ.0) TRIID = FACES_TAG(I,5)

	!	x1 = NODES(FACES_TAG(I,1),1)
	!	x2 = NODES(FACES_TAG(I,2),1)
	!	y1 = NODES(FACES_TAG(I,1),2)
	!	y2 = NODES(FACES_TAG(I,2),2)

	!	IF (x1.GT.x2) THEN
	!		x3 = x1
	!		x1 = x2
	!		x2 = x3
	!	END IF

	!	DO J  = 1, NP_FC
	!		yp = 0.5d0*(y1+y2)
	!		xp = 0.5d0*FACES_INF(I)*(QUA_ABS_LINE(J)+1.d0) + x1

	!		x(1) = NODES(TRIANGLES_TAG(TRIID,1),1)
	!		y(1) = NODES(TRIANGLES_TAG(TRIID,1),2)
	!		x(2) = NODES(TRIANGLES_TAG(TRIID,2),1)
	!		y(2) = NODES(TRIANGLES_TAG(TRIID,2),2)
	!		x(3) = NODES(TRIANGLES_TAG(TRIID,3),1)
	!		y(3) = NODES(TRIANGLES_TAG(TRIID,3),2)
	!		x(4) = x(1)
	!		y(4) = y(1)

	!		A=(x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1))
	!		B=y(3)-y(1)
	!		C=x(1)-x(3)
	!		D=(x(3)-x(1))*y(1)-(y(3)-y(1))*x(1)

	!		E=(x(3)-x(1))*(y(2)-y(1))-(x(2)-x(1))*(y(3)-y(1))
	!		F=y(2)-y(1)
	!		G=x(1)-x(2)
	!		H=(x(2)-x(1))*y(1)-(y(2)-y(1))*x(1)

	!		xi = B/A*xp+C/A*yp+D/A
	!		eta = F/E*xp+G/E*yp+H/E

	!		pm = 0.d0
	!		DO M = 1, NDOF_TRI
	!			L = 1
	!			DO KL = 0,DEG
	!				DO JL = 0, KL
	!					IL = KL - JL
	!					pm(M) = pm(M) + NODFUN_TRI_REF(M,L)*(xi**REAL(IL,DBL))*(eta**REAL(JL,DBL))
	!					L = L + 1
	!				END DO
	!			END DO
	!		END DO
	!		NODFUN_QUA_P_LINE(J,:,I) = pm
	!	END DO
	!END IF
	!END DO
	!END IF

	DEALLOCATE (WEIG,ABSC,pm,pl)

	END SUBROUTINE Init_Integration

END MODULE
