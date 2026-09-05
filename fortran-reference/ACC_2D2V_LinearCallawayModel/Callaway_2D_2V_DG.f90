PROGRAM Callaway_2D_2V_DG
! Solving the 2D2V TIME INDEPEDENT Callaway model equation using DG method
! no large linear system
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
IMPLICIT NONE
INTEGER :: ios, STEP
INTEGER*8 :: t1,t2,count_rate,count_max
REAL (KIND=DBL) :: RESIDUL,MASS0
CHARACTER(LEN=50) :: FILENAME
CHARACTER(LEN=4) :: FLAG

CALL MKL_SET_NUM_THREADS (1)
CALL MKL_SET_DYNAMIC (0)
CALL OMP_SET_NUM_THREADS (24)

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

!MASS0=MASS
WRITE (*,*) '**** Init End ****'

CALL SYSTEM_CLOCK (t1,count_rate,count_max)

STEP = 1
DO
	!*******************************************
	!DVM
	!CALL Calculate_FLUX_WALL()
	CALL DG_Solver_VDF ()
	IF (ACCFLAG.EQ.0)	CALL Calculate_Macro_Properties ()

	!*******************************************
	!Acceleratioin
	IF (ACCFLAG.EQ.1) THEN
		CALL Calculate_SRC_ACC_HoTfromDVM ()
		CALL Global_Problem_Solver_ACC ()
		CALL Local_Problem_Solver_ACC ()
		CALL Correct_VDF_Calculate_Macro_Properties ()
	END IF

    CALL Calculate_Residual_T (RESIDUL)
	WRITE (*,*) STEP, RESIDUL, MASS

	IF (RESIDUL.LT.TOL) EXIT
	IF (STEP.GE.TMAX) EXIT
	STEP = STEP + 1
END DO

IF (ACCFLAG.EQ.1) CALL Release_PARDISO ()
CALL SYSTEM_CLOCK (t2,count_rate,count_max)

!IF (ACCFLAG.EQ.0) CALL Compress_VDF()
CALL Out_Put_Conduction ()

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
END PROGRAM Callaway_2D_2V_DG
