MODULE SOLVER
USE CONSTANT
USE GLOBAL_VARIABLE
USE MATRICES
USE VELOCITY_DISTRIBUTION_FUNCTION
USE VELOCITY_GRID
USE SPATIAL_GRID
USE BOUNDARY_CONDITION
USE INTEGRATIONS
IMPLICIT NONE

CONTAINS

SUBROUTINE DG_Solver_VDF ()
IMPLICIT NONE
EXTERNAL:: DGETRF, DGETRS
INTEGER :: J1,J2,I,M,IL,L
INTEGER :: TRID, TRIDext, FCID
INTEGER :: BCID, JS1, JS2
INTEGER :: PAIR,INFO
REAL (KIND = DBL) :: speed


!$OMP PARALLEL DO PRIVATE (J1,J2,I,TRID,TRIDext,M,L,IL,INFO,BCID,FCID,JS1,JS2, &
!$OMP PAIR,speed,IPIV,A_SOL,A_SRC,FW)
DO J2 = 1, NAZIM
	DO J1 = 1, NPOLE
!**************************************************************
		DO I = 1, N_TRIS
			TRID = TRI_ORDER(I,J1,J2)
			!calculate A_SOL
			A_SOL = 0.d0
			DO L = 1, NDOF_TRI
    			DO M = 1, NDOF_TRI
    				A_SOL(M,L) = 1.d0/TAU_C*INT_NODFUNC_TRI_TRI(M,L,TRID) &
    						- CX(J1,J2)*INT_NODFUNC_TRI_TRI_X(TRID,M,L)-CY(J1,J2)*INT_NODFUNC_TRI_TRI_Y(TRID,M,L)
    				DO IL = 1,3
    					speed = CX(J1,J2)*TRIANGLES_INF(TRID,2*IL)+CY(J1,J2)*TRIANGLES_INF(TRID,2*IL+1)
    					A_SOL(M,L) = A_SOL(M,L) + 0.5d0*(speed+ABS(speed))*INT_NODFUNC_TRI_FC(TRID,IL,M,L)
    				END DO
    			END DO
    		END DO

!*****************************************************************************************************
			!calculate A_SRC
			A_SRC = 0.d0
			DO M = 1, NDOF_TRI
				DO L = 1, NDOF_TRI
					A_SRC(M) = A_SRC(M) + Cv*T_s(L,TRID)/4.d0/PI/TAU_R*INT_NODFUNC_TRI_TRI(L,M,TRID) &
						+ 1.d0/TAU_N*INT_NODFUNC_TRI_TRI(L,M,TRID) &
						* (Cv*T_s(L,TRID)/4.d0/PI+3.d0/4.d0/PI*(Qx_s(L,TRID)*CX(J1,J2)+Qy_s(L,TRID)*CY(J1,J2))/Vg/Vg)
				END DO
			END DO

			DO IL = 1, 3
				speed = CX(J1,J2)*TRIANGLES_INF(TRID,2*IL)+CY(J1,J2)*TRIANGLES_INF(TRID,2*IL+1)
				FCID = TRIANGLES_TAG(TRID,3+IL)
				BCID = FACES_TAG(FCID,3)

				IF (BCID.EQ.0) THEN
					!interior
					IF (FACES_TAG(FCID,5).EQ.TRID) TRIDext = FACES_TAG(FCID,6)  !neighbor triangle
					IF (FACES_TAG(FCID,6).EQ.TRID) TRIDext = FACES_TAG(FCID,5)

					DO L = 1, NDOF_TRI
						A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*VDF(L,TRIDext,J1,J2)
					END DO
				ELSE
					!boundary
					!IF (BC_TYP(BCID).EQ.1) THEN
						!thermalisation
						FW = Cv/4.d0/PI*BC_TEMP(BCID)*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)
						A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*FW
					!END IF
					!IF (BC_TYP(BCID).EQ.2) THEN
						!non-thermalisation
					!	FW = -FLUX_WALL(:,FCID)
					!	A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*FW
					!END IF

					!IF (BC_TYP(BCID).EQ.3) THEN
						!Period
					!	PAIR = FACES_TAG(FCID,4)
					!	IF (FACES_TAG(PAIR,5).EQ.0) TRIDext = FACES_TAG(PAIR,6)
					!	IF (FACES_TAG(PAIR,6).EQ.0) TRIDext = FACES_TAG(PAIR,5)

					!	DO L = 1, NDOF_TRI
					!		A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*VDF(L,TRIDext,J1,J2)
					!	END DO
					!END IF
				END IF
			END DO
!**************************************************************************************************

			!local solver
			CALL DGETRF (NDOF_TRI,NDOF_TRI,A_SOL,NDOF_TRI,IPIV,INFO)
    		!IF (INFO.NE.0) WRITE (*,*) 'Warning LU factorization error: ',INFO
    		CALL DGETRS ('N',NDOF_TRI,1,A_SOL,NDOF_TRI,IPIV,A_SRC,NDOF_TRI,INFO)

			VDF(:,TRID,J1,J2) = A_SRC
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!STOP

END SUBROUTINE DG_Solver_VDF

END MODULE
