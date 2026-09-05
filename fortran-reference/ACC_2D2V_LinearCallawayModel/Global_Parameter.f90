!--------------------------------------------

MODULE CONSTANT
!Define global constants
IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND (p=5)
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND (p=15)
INTEGER, PARAMETER :: INTL=SELECTED_INT_KIND (12)

REAL (KIND=DBL), PARAMETER :: PI=4.d0*atan(1.d0)        !Pi
COMPLEX (KIND=DBL), PARAMETER :: IMG_I = CMPLX(0.0,1.0,DBL) !imaginary unit

REAL (KIND=DBL), PARAMETER :: BOLTZ=1.3806485279d-23    !Boltzmann Constant
REAL (KIND=DBL), PARAMETER :: AFGDL=6.022140857d+23     !Avogadro constant

REAL (KIND=DBL), PARAMETER :: MIN_ERR=1.0d-16
REAL (KIND=DBL), PARAMETER :: MID_ERR=1.0d-10
REAL (KIND=DBL), PARAMETER :: MAX_ERR=1.0d-6

END MODULE CONSTANT

!-----------------------------------------------

MODULE GLOBAL_VARIABLE
!Define and initial global variables
USE CONSTANT
IMPLICIT NONE

REAL (KIND=DBL), SAVE :: TOL=1.0d-5                            !iteration tolerance
INTEGER, SAVE :: TMAX=10000                                 !maximum iteration steps
INTEGER, SAVE :: NPOLE=20, NAZIM=12                         !number of discretized polar and azimuthal angles
INTEGER, SAVE :: DEG=1                                   !degree of polynomial of DG discretization
INTEGER, SAVE :: NDOF_TRI                               !number of degree of freedom in element (triangle): NDEG_TRI=(DEG+1)*(DEG+2)/2
INTEGER, SAVE :: NDOF_FC                               !number of degree of freedom on face: NDEG_FC=DEG+1
INTEGER, SAVE :: NP_TRI                                !number of quadrature points in element
INTEGER, SAVE :: NP_FC                                 !number of quadrature points on face

REAL (KIND=DBL), SAVE :: Cv=0.d0
REAL (KIND=DBL), SAVE :: Vg=1.d0
REAL (KIND=DBL), SAVE :: TAU_R=0.0
REAL (KIND=DBL), SAVE :: TAU_N=0.0
REAL (KIND=DBL), SAVE :: TAU_C=0.0
REAL (KIND=DBL), SAVE :: TAU_THR=1.d0

CHARACTER (LEN=50), SAVE :: FNAME_MSH='Gmsh.msh'                   !File name of spatial mesh

INTEGER, SAVE :: NBC = 1                                   !# of boundary condition
INTEGER, SAVE :: ACCFLAG                     !0-no acceleration, 1-acceleration

CONTAINS
   SUBROUTINE Init_Global_Variables (fin)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: fin
   INTEGER :: ios

   NAMELIST /ITERATION/TOL,TMAX
   NAMELIST /VELMSH/NPOLE,NAZIM
   NAMELIST /DG/DEG
   NAMELIST /FLOW/Cv,Vg,TAU_R,TAU_N,TAU_THR
   NAMELIST /FILENAME/FNAME_MSH
   NAMELIST /N_BC/NBC
   NAMELIST /GSIS/ACCFLAG

   WRITE (*,*) '**** Initial global variables ****'
   READ (fin, NML=ITERATION, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read ITERATION parameters: ', ios
   WRITE (*,100) 'Iteration Tolerance: ', TOL
   100 FORMAT (1X, A30, ES10.3)
   WRITE (*,101) 'Maximum iteration step: ', TMAX
   101 FORMAT (1X, A30, I10)
   WRITE (*,*)

   READ (fin, NML=GSIS, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read GSIS parameters: ', ios
   IF (ACCFLAG.EQ.0) WRITE (*,105) 'Iteration Scheme: ','CIS'
   IF (ACCFLAG.EQ.1) WRITE (*,105) 'Iteration Scheme: ','GSIS'
   105 FORMAT (1X, A30, A4)
   WRITE (*,*)

   READ (fin, NML=VELMSH, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read VELMESH parameters: ', ios
   NAZIM = (NAZIM/2)*2
   WRITE (*,102) '# of polar angles: ', NPOLE
   102 FORMAT (1X,A30,I5)
   WRITE (*,102) '# of azimuthal angles: ', NAZIM
   WRITE (*,*)

   READ (fin, NML=DG, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read DG parameters: ', ios
   NDOF_TRI=(DEG+1)*(DEG+2)/2
   NDOF_FC=DEG+1

   IF(DEG.EQ.1) NP_TRI = 3
   IF(DEG.EQ.2) NP_TRI = 6
   IF(DEG.EQ.3) NP_TRI = 12
   IF(DEG.EQ.4) NP_TRI = 16

   NP_FC=15
   WRITE (*,103) 'Degree of polynomial: ', DEG
   103 FORMAT (1X,A30,I5)
   WRITE (*,103) 'Degree of freedom in element: ', NDOF_TRI
   WRITE (*,103) 'Degree of freedom on face: ', NDOF_FC
   WRITE (*,103) 'Quadrature points in element: ', NP_TRI
   WRITE (*,103) 'Quadrature points on face: ', NP_FC
   WRITE (*,*)

   READ (fin, NML=FLOW, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read FLOW parameters: ', ios
   WRITE (*,104)'Volumetric specific heat: ', Cv
   104 FORMAT (1X, A30, F10.3)
   WRITE (*,104) 'Magnitude of group velocity: ', Vg

   TAU_C = 1.d0/(1.d0/TAU_R+1.d0/TAU_N)
   WRITE (*,104)'tau of resistive scattering: ', TAU_R
   WRITE (*,104)'tau of normal scattering: ', TAU_N
   WRITE (*,104)'conbinating tau: ', TAU_C
   WRITE (*,*)


   READ (fin, NML=FILENAME, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read FILENAME parameters: ', ios
   WRITE (*,109) 'File for spatial mesh: ', FNAME_MSH
   109 FORMAT (1X,A30,A20)
   WRITE (*,*)

   READ (fin, NML=N_BC, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read number of boundary condtion: ', ios
   WRITE (*,102) '# of boundary conditions: ', NBC
   WRITE (*,*)


   END SUBROUTINE Init_Global_Variables

END MODULE GLOBAL_VARIABLE
