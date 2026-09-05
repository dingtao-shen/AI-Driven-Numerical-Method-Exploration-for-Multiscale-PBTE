PROGRAM Callaway_2D_2V_DG
! Solving the 2D2V TIME INDEPEDENT Callaway model equation using DG method
! no large linear system
!
! ---------------------------------------------------------------------------
! INSTRUMENTED COPY of fortran-reference/.../Callaway_2D_2V_DG.f90
!
! The numerics are byte-for-byte the reference driver.  What is added:
!   * CALL DUMP_INIT / DUMP_* stage snapshots  (no-ops unless PYBTE_DUMP_DIR
!     is set in the environment)
!   * residual_history.txt, so Gate 0 has a machine-readable golden log
! What is *changed*:
!   * OMP_SET_NUM_THREADS(24) -> OMP_SET_NUM_THREADS(1); the reference is
!     built serially for bit-reproducibility (Proposal 1 §Phase 0.5).
! Nothing else.
! ---------------------------------------------------------------------------
USE CONSTANT
USE GLOBAL_VARIABLE
USE BOUNDARY_CONDITION
USE SPATIAL_GRID
USE VELOCITY_GRID
USE NODAL_FUNCTION
USE INTEGRATIONS
USE VELOCITY_DISTRIBUTION_FUNCTION
USE MATRICES
USE ACCELERATION
USE SOLVER
USE OUT_PUT
USE DUMP_UTILS
IMPLICIT NONE
EXTERNAL :: DUMP_R8, DUMP_I4
INTEGER :: ios, STEP
INTEGER*8 :: t1,t2,count_rate,count_max
REAL (KIND=DBL) :: RESIDUL
CHARACTER(LEN=50) :: FILENAME
CHARACTER(LEN=4) :: FLAG
LOGICAL :: WANT

CALL MKL_SET_NUM_THREADS (1)
CALL MKL_SET_DYNAMIC (0)
CALL OMP_SET_NUM_THREADS (1)

CALL DUMP_INIT ()

!Initial Global Variables and boundary condition
OPEN (UNIT=10, FILE='./control.in', STATUS='OLD', ACTION='READ', IOSTAT=ios)
IF (ios.NE.0) THEN
   	WRITE (*,*) 'Failed to open control file: ', ios
ELSE
	CALL Init_Global_Variables (10)
	CALL Init_Boundary_Conditions (10)
END IF
CLOSE (UNIT=10)


CALL Init_Spatial_Grid ()    !Initial spatial mesh
CALL Init_Velocity_Grid ()   !Initial velocity mesh
CALL Init_Triangle_Order ()  !Initial Triangle order
CALL Init_Basis_Function ()  !Initial Basis
CALL Init_Integration ()     !Initial integration of basis functions
CALL Init_Velocity_Distribution_Function ()    !Initial velocity distribution functions
!CALL Calculate_Macro_Properties ()
!************************************************
IF (ACCFLAG.EQ.1) THEN
	CALL Init_Acceleration_New ()
	CALL Init_PARDISO ()
END IF
!************************************************

IF (DUMP_ON) CALL DUMP_SETUP ()

!MASS0=MASS
WRITE (*,*) '**** Init End ****'

OPEN (15,FILE='./residual_history.txt',STATUS='REPLACE',ACTION='WRITE')
WRITE (15,'(A)') '# step residual mass'

CALL SYSTEM_CLOCK (t1,count_rate,count_max)

STEP = 1
DO
	WANT = DUMP_WANT_ITER (STEP)
	!*******************************************
	!DVM
	!CALL Calculate_FLUX_WALL()
	CALL DG_Solver_VDF ()
	IF (WANT) CALL DUMP_STAGE_VDF (STEP, 'vdf_after_sweep')
	IF (ACCFLAG.EQ.0)	CALL Calculate_Macro_Properties ()
	IF (WANT.AND.(ACCFLAG.EQ.0)) CALL DUMP_STAGE_MOMENTS (STEP, 'after_moments')

	!*******************************************
	!Acceleratioin
	IF (ACCFLAG.EQ.1) THEN
		CALL Calculate_SRC_ACC_HoTfromDVM ()
		IF (WANT) CALL DUMP_STAGE_AASRC (STEP)
		CALL Global_Problem_Solver_ACC ()
		IF (WANT) CALL DUMP_STAGE_TRACE (STEP)
		CALL Local_Problem_Solver_ACC ()
		IF (WANT) CALL DUMP_STAGE_UQ (STEP)
		CALL Correct_VDF_Calculate_Macro_Properties ()
		IF (WANT) THEN
			CALL DUMP_STAGE_VDF (STEP, 'vdf_after_correction')
			CALL DUMP_STAGE_MOMENTS (STEP, 'after_correction')
		END IF
	END IF

    CALL Calculate_Residual_T (RESIDUL)
	WRITE (*,*) STEP, RESIDUL, MASS
	WRITE (15,'(I10,2ES26.16E3)') STEP, RESIDUL, MASS

	IF (RESIDUL.LT.TOL) EXIT
	IF (STEP.GE.TMAX) EXIT
	STEP = STEP + 1
END DO

CLOSE (15)

IF (ACCFLAG.EQ.1) CALL Release_PARDISO ()
CALL SYSTEM_CLOCK (t2,count_rate,count_max)

!IF (ACCFLAG.EQ.0) CALL Compress_VDF()
CALL Out_Put_Conduction ()

IF (DUMP_ON) CALL DUMP_FINAL (STEP)

WRITE (*,100) REAL(t2-t1)/REAL(count_rate)
100 FORMAT (1X,'**** RUN END, Time = ', F10.3, ' [s] ****')

IF (ACCFLAG.EQ.0) FLAG=' CIS'
IF (ACCFLAG.EQ.1) FLAG='GSIS'
FILENAME = './RunTime.txt'
OPEN (14,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
WRITE(14,101)FLAG, TAU_R, TAU_N, N_TRIS, NPOLE, NAZIM, Step, REAL(t2-t1)/REAL(count_rate)
CLOSE(14)

101 FORMAT (1X,A4,2ES10.2,3I6,I8,F20.1)

STOP

CONTAINS

!==========================================================================
SUBROUTINE DUMP_SETUP ()
IMPLICIT NONE
INTEGER :: D(4)

! ---- scalars ----------------------------------------------------------
CALL DUMP_SCALAR_I4 ('s_N_NDS',   N_NDS)
CALL DUMP_SCALAR_I4 ('s_N_ELES',  N_ELES)
CALL DUMP_SCALAR_I4 ('s_N_TRIS',  N_TRIS)
CALL DUMP_SCALAR_I4 ('s_N_FCS',   N_FCS)
CALL DUMP_SCALAR_I4 ('s_N_FCS_B', N_FCS_B)
CALL DUMP_SCALAR_I4 ('s_N_FCS_I', N_FCS_I)
CALL DUMP_SCALAR_I4 ('s_N_FCS_WALL_SYM', N_FCS_WALL_SYM)
CALL DUMP_SCALAR_I4 ('s_NDOF_TRI', NDOF_TRI)
CALL DUMP_SCALAR_I4 ('s_NDOF_FC',  NDOF_FC)
CALL DUMP_SCALAR_I4 ('s_NP_TRI',   NP_TRI)
CALL DUMP_SCALAR_I4 ('s_NP_FC',    NP_FC)
CALL DUMP_SCALAR_I4 ('s_NPOLE',    NPOLE)
CALL DUMP_SCALAR_I4 ('s_NAZIM',    NAZIM)
CALL DUMP_SCALAR_I4 ('s_DEG',      DEG)
CALL DUMP_SCALAR_I4 ('s_ACCFLAG',  ACCFLAG)
CALL DUMP_SCALAR_R8 ('s_TAU_R',    TAU_R)
CALL DUMP_SCALAR_R8 ('s_TAU_N',    TAU_N)
CALL DUMP_SCALAR_R8 ('s_TAU_C',    TAU_C)
CALL DUMP_SCALAR_R8 ('s_TAU_THR',  TAU_THR)
CALL DUMP_SCALAR_R8 ('s_Cv',       Cv)
CALL DUMP_SCALAR_R8 ('s_Vg',       Vg)
CALL DUMP_SCALAR_R8 ('s_Hmin',     Hmin)
CALL DUMP_SCALAR_R8 ('s_TOL',      TOL)

! ---- velocity mesh ----------------------------------------------------
D(1) = NPOLE
CALL DUMP_R8 ('THE',  THE,  NPOLE, 1, D)
CALL DUMP_R8 ('WTHE', WTHE, NPOLE, 1, D)
D(1) = NAZIM
CALL DUMP_R8 ('PHI',  PHI,  NAZIM, 1, D)
CALL DUMP_R8 ('WPHI', WPHI, NAZIM, 1, D)
D(1) = NPOLE
D(2) = NAZIM
CALL DUMP_R8 ('DOMEGA', DOMEGA, NPOLE*NAZIM, 2, D)
CALL DUMP_R8 ('CX',     CX,     NPOLE*NAZIM, 2, D)
CALL DUMP_R8 ('CY',     CY,     NPOLE*NAZIM, 2, D)

! ---- spatial mesh -----------------------------------------------------
D(1) = N_NDS
D(2) = 2
CALL DUMP_R8 ('NODES', NODES, N_NDS*2, 2, D)
D(1) = N_TRIS
D(2) = 6
CALL DUMP_I4 ('TRIANGLES_TAG', TRIANGLES_TAG, N_TRIS*6, 2, D)
D(2) = 7
CALL DUMP_R8 ('TRIANGLES_INF', TRIANGLES_INF, N_TRIS*7, 2, D)
D(1) = N_TRIS
D(2) = 1
CALL DUMP_R8 ('TRIANGLES_Hmin', TRIANGLES_Hmin, N_TRIS, 1, D)
D(1) = N_FCS
D(2) = 8
CALL DUMP_I4 ('FACES_TAG', FACES_TAG, N_FCS*8, 2, D)
D(1) = N_FCS
CALL DUMP_R8 ('FACES_INF', FACES_INF, N_FCS, 1, D)

! ---- sweep ordering ---------------------------------------------------
D(1) = N_TRIS
D(2) = NPOLE
D(3) = NAZIM
CALL DUMP_I4 ('TRI_ORDER', TRI_ORDER, N_TRIS*NPOLE*NAZIM, 3, D)

! ---- basis ------------------------------------------------------------
D(1) = NDOF_FC
D(2) = NDOF_FC
CALL DUMP_R8 ('NODFUN_FC_REF', NODFUN_FC_REF, NDOF_FC*NDOF_FC, 2, D)
D(1) = NDOF_FC
CALL DUMP_R8 ('NODE_FC_REF', NODE_FC_REF, NDOF_FC, 1, D)
D(1) = NDOF_TRI
D(2) = NDOF_TRI
CALL DUMP_R8 ('NODFUN_TRI_REF', NODFUN_TRI_REF, NDOF_TRI*NDOF_TRI, 2, D)
D(1) = NDOF_TRI
D(2) = 2
CALL DUMP_R8 ('NODE_TRI_REF', NODE_TRI_REF, NDOF_TRI*2, 2, D)

! ---- integration tensors ---------------------------------------------
D(1) = NDOF_TRI
D(2) = N_TRIS
CALL DUMP_R8 ('INT_NODFUNC_TRI', INT_NODFUNC_TRI, NDOF_TRI*N_TRIS, 2, D)

D(1) = NDOF_TRI
D(2) = NDOF_TRI
D(3) = N_TRIS
CALL DUMP_R8 ('INT_NODFUNC_TRI_TRI', INT_NODFUNC_TRI_TRI, &
              NDOF_TRI*NDOF_TRI*N_TRIS, 3, D)

D(1) = NDOF_TRI
D(2) = NDOF_TRI
D(3) = NDOF_TRI
D(4) = N_TRIS
CALL DUMP_R8 ('INT_NODFUNC_TRI_TRI_TRI', INT_NODFUNC_TRI_TRI_TRI, &
              NDOF_TRI*NDOF_TRI*NDOF_TRI*N_TRIS, 4, D)

D(1) = N_TRIS
D(2) = NDOF_TRI
D(3) = NDOF_TRI
CALL DUMP_R8 ('INT_NODFUNC_TRI_TRI_X', INT_NODFUNC_TRI_TRI_X, &
              N_TRIS*NDOF_TRI*NDOF_TRI, 3, D)
CALL DUMP_R8 ('INT_NODFUNC_TRI_TRI_Y', INT_NODFUNC_TRI_TRI_Y, &
              N_TRIS*NDOF_TRI*NDOF_TRI, 3, D)

D(1) = N_TRIS
D(2) = 3
D(3) = NDOF_TRI
D(4) = NDOF_TRI
CALL DUMP_R8 ('INT_NODFUNC_TRI_FC', INT_NODFUNC_TRI_FC, &
              N_TRIS*3*NDOF_TRI*NDOF_TRI, 4, D)
CALL DUMP_R8 ('INT_NODFUNC_TRI_TRI_FC', INT_NODFUNC_TRI_TRI_FC, &
              N_TRIS*3*NDOF_TRI*NDOF_TRI, 4, D)

D(1) = NDOF_TRI
D(2) = NDOF_FC
D(3) = N_TRIS
D(4) = 3
CALL DUMP_R8 ('INT_NODFUNC_TRI_FC_FC', INT_NODFUNC_TRI_FC_FC, &
              NDOF_TRI*NDOF_FC*N_TRIS*3, 4, D)

D(1) = NDOF_FC
D(2) = NDOF_FC
D(3) = N_FCS
CALL DUMP_R8 ('INT_NODFUNC_FC_FC', INT_NODFUNC_FC_FC, &
              NDOF_FC*NDOF_FC*N_FCS, 3, D)

D(1) = NDOF_FC
D(2) = N_FCS
CALL DUMP_R8 ('INT_NODFUNC_FC', INT_NODFUNC_FC, NDOF_FC*N_FCS, 2, D)

! ---- element quadrature ----------------------------------------------
D(1) = NP_TRI
CALL DUMP_R8 ('QUA_ABS_X', QUA_ABS_X, NP_TRI, 1, D)
CALL DUMP_R8 ('QUA_ABS_Y', QUA_ABS_Y, NP_TRI, 1, D)
CALL DUMP_R8 ('QUA_WEI',   QUA_WEI,   NP_TRI, 1, D)
D(1) = NDOF_TRI
D(2) = NP_TRI
CALL DUMP_R8 ('NODFUN_QUA_P', NODFUN_QUA_P, NDOF_TRI*NP_TRI, 2, D)

! ---- GSIS setup -------------------------------------------------------
IF (ACCFLAG.EQ.1) THEN
   D(1) = 3
   CALL DUMP_R8 ('ST', ST, 3, 1, D)

   D(1) = 7*NDOF_TRI
   D(2) = 7*NDOF_TRI
   D(3) = N_TRIS
   CALL DUMP_R8 ('inv_AA_SOL', inv_AA_SOL, 7*NDOF_TRI*7*NDOF_TRI*N_TRIS, 3, D)

   D(1) = 7*NDOF_TRI
   D(2) = 3*NDOF_FC
   D(3) = 3
   D(4) = N_TRIS
   CALL DUMP_R8 ('AA_TRACE', AA_TRACE, 7*NDOF_TRI*3*NDOF_FC*3*N_TRIS, 4, D)

   D(1) = 3*NDOF_FC
   D(2) = 7*NDOF_TRI
   D(3) = 3
   D(4) = N_TRIS
   CALL DUMP_R8 ('BA_SOL', BA_SOL, 3*NDOF_FC*7*NDOF_TRI*3*N_TRIS, 4, D)

   D(1) = SIZE(KKAcomp)
   CALL DUMP_R8 ('KKAcomp', KKAcomp, SIZE(KKAcomp), 1, D)
   D(1) = SIZE(JKKA)
   CALL DUMP_I4 ('JKKA', JKKA, SIZE(JKKA), 1, D)
   D(1) = SIZE(IKKA)
   CALL DUMP_I4 ('IKKA', IKKA, SIZE(IKKA), 1, D)
END IF

END SUBROUTINE DUMP_SETUP

!==========================================================================
SUBROUTINE DUMP_STAGE_VDF (S, TAG)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
CHARACTER (LEN=*), INTENT(IN) :: TAG
INTEGER :: D(4)
IF (.NOT.DUMP_VDF_ON) RETURN
D(1) = NDOF_TRI
D(2) = N_TRIS
D(3) = NPOLE
D(4) = NAZIM
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,TAG)), VDF, NDOF_TRI*N_TRIS*NPOLE*NAZIM, 4, D)
END SUBROUTINE DUMP_STAGE_VDF

!==========================================================================
SUBROUTINE DUMP_STAGE_MOMENTS (S, TAG)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
CHARACTER (LEN=*), INTENT(IN) :: TAG
INTEGER :: D(2)
D(1) = NDOF_TRI
D(2) = N_TRIS
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Ts_'//TAG)),  T_s,  NDOF_TRI*N_TRIS, 2, D)
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Qxs_'//TAG)), Qx_s, NDOF_TRI*N_TRIS, 2, D)
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Qys_'//TAG)), Qy_s, NDOF_TRI*N_TRIS, 2, D)
D(1) = N_TRIS
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Temp_'//TAG)), Temp, N_TRIS, 1, D)
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Qx_'//TAG)),   Qx,   N_TRIS, 1, D)
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'Qy_'//TAG)),   Qy,   N_TRIS, 1, D)
END SUBROUTINE DUMP_STAGE_MOMENTS

!==========================================================================
SUBROUTINE DUMP_STAGE_AASRC (S)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
INTEGER :: D(2)
D(1) = 7*NDOF_TRI
D(2) = N_TRIS
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'AA_SRC')), AA_SRC, 7*NDOF_TRI*N_TRIS, 2, D)
END SUBROUTINE DUMP_STAGE_AASRC

!==========================================================================
SUBROUTINE DUMP_STAGE_TRACE (S)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
INTEGER :: D(2)
D(1) = 3*NDOF_FC
D(2) = N_FCS
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'U_TRACE')), U_TRACE, 3*NDOF_FC*N_FCS, 2, D)
D(1) = N_FCS*3*NDOF_FC
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'FFA')), FFA, N_FCS*3*NDOF_FC, 1, D)
END SUBROUTINE DUMP_STAGE_TRACE

!==========================================================================
SUBROUTINE DUMP_STAGE_UQ (S)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
INTEGER :: D(2)
D(1) = 7*NDOF_TRI
D(2) = N_TRIS
CALL DUMP_R8 (TRIM(DUMP_ITNAME(S,'UQ')), UQ, 7*NDOF_TRI*N_TRIS, 2, D)
END SUBROUTINE DUMP_STAGE_UQ

!==========================================================================
SUBROUTINE DUMP_FINAL (S)
IMPLICIT NONE
INTEGER, INTENT(IN) :: S
INTEGER :: D(4)
CALL DUMP_SCALAR_I4 ('final_steps', S)
D(1) = N_TRIS
CALL DUMP_R8 ('final_Temp', Temp, N_TRIS, 1, D)
CALL DUMP_R8 ('final_Qx',   Qx,   N_TRIS, 1, D)
CALL DUMP_R8 ('final_Qy',   Qy,   N_TRIS, 1, D)
D(1) = NDOF_TRI
D(2) = N_TRIS
CALL DUMP_R8 ('final_Ts',  T_s,  NDOF_TRI*N_TRIS, 2, D)
CALL DUMP_R8 ('final_Qxs', Qx_s, NDOF_TRI*N_TRIS, 2, D)
CALL DUMP_R8 ('final_Qys', Qy_s, NDOF_TRI*N_TRIS, 2, D)
IF (DUMP_VDF_ON) THEN
   D(1) = NDOF_TRI
   D(2) = N_TRIS
   D(3) = NPOLE
   D(4) = NAZIM
   CALL DUMP_R8 ('final_VDF', VDF, NDOF_TRI*N_TRIS*NPOLE*NAZIM, 4, D)
END IF
END SUBROUTINE DUMP_FINAL

END PROGRAM Callaway_2D_2V_DG
