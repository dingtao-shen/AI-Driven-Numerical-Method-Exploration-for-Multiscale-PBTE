MODULE BOUNDARY_CONDITION
!Define boundary condition
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE

CHARACTER (LEN=20), ALLOCATABLE, DIMENSION(:), SAVE :: BC_NAME     !NAME
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: BC_PHYID               !Physical id associated with the mesh
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: BC_TYP             !Boundary type: 1-wall, 2-pressure inlet/outlet, 3-period, 4-symmetry,
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: BC_TEMP !perturbed temperature
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: BC_XOFF, BC_YOFF  !offset for periodic condition
INTEGER, SAVE :: BCID_MASTER = 0, BCID_SLAVE = 0

CONTAINS

	SUBROUTINE Init_Boundary_Conditions (fin)
	IMPLICIT NONE
	INTEGER, INTENT (IN) :: fin
	INTEGER :: ios, I

	NAMELIST /BC/BC_NAME,BC_PHYID,BC_TYP,BC_TEMP,BC_XOFF,BC_YOFF

	WRITE (*,*) '**** Initial BC parameters ****'
	ALLOCATE (BC_NAME(NBC),BC_PHYID(NBC),BC_TYP(NBC),BC_TEMP(NBC),BC_XOFF(NBC),BC_YOFF(NBC))

	READ (fin, NML=BC, IOSTAT=ios)
	IF (ios.NE.0) WRITE(*,*) 'Error to read BC parameters: ', ios

	WRITE (*,100) 'BC_ID','NAME','PHY_ID','TYPE','TEMP','XOFF','YOFF'
	100 FORMAT (1X,A7,A20,2A7,3A11)
	DO I = 1, NBC
		WRITE (*,101) I,TRIM(BC_NAME(I)),BC_PHYID(I),BC_TYP(I),BC_TEMP(I),BC_XOFF(I),BC_YOFF(I)
		101 FORMAT (1X,I7,A20,2I7,3ES11.3)
	END DO
	WRITE (*,*)

	!FIND BC index of the periodic condition
	DO I = 1, NBC
		IF (TRIM(BC_NAME(I))=='Slave') BCID_SLAVE = I
		IF (TRIM(BC_NAME(I))=='Master') BCID_MASTER = I
	END DO
	WRITE (*,107) BCID_MASTER, BCID_SLAVE
	107 FORMAT (1X,'The BC index of periodic master and slave are: ',2I5)
	WRITE (*,*)

	END SUBROUTINE Init_Boundary_Conditions

END MODULE BOUNDARY_CONDITION
