MODULE SPATIAL_GRID
!Define spatial mesh
USE CONSTANT
USE GLOBAL_VARIABLE
USE BOUNDARY_CONDITION
IMPLICIT NONE

INTEGER, SAVE :: N_NDS, N_ELES !# of nodes and elements in gmsh file
INTEGER, SAVE :: N_FCS, N_FCS_B, N_FCS_I,N_TRIS !# of faces and triangles in the mesh
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: NODES !dimension (N_NDS,2): (:,1)/(:,2) are x/y of i-th node
INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE :: TRIANGLES_TAG ! dimension (N_TRIS,6): (:,1:3)-node index
                                                            ! 3
                                                            ! | \        edge 1: 1--2
                                                            ! |  \		 edge 2: 2--3
                                                            ! 1---2      edge 3: 3--1
                                                            ! (:,4:6)-faces index for each edge
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: TRIANGLES_INF ! dimension (N_TRIS,7): (:,1)-area
                                                                    ! (:,2:3)-outward normal vector of edge 1
                                                                    ! (:,4:5)-outward normal vector of edge 2
                                                                    ! (:,6:7)-outward normal vector of edge 3
INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE :: FACES_TAG ! dimension (N_FCS, 8): (:,1:2)-node index
                                                        ! (:,3)-BCflage, >0 is the index of BC condition, =0 is the interior faces
                                                        ! (:,4)-index of periodic pair if it has one,
                                                        !       or the index of wall and symmetric faces in the wall and symmetric queue
                                                        ! (:, 5)-index of triangle(+), (:,6)-index of triangle(-)
                                                        ! (:, 7)-local face ID in triangle(+)
                                                        ! (:,8)-local face ID in triangle(-)
                                                        !                2
                                                        !                |
                                                        !   triangle(+)  |  triangle(-)
                                                        !	             |
                                                        !                1
                                                        !curl right hand fingers from 1 to 2, right thumb points to the triangle(-)

REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: FACES_INF   !dimension (N_FCS): length of faces
INTEGER, SAVE :: N_FCS_WALL_SYM  ! # of wall boundary faces, and symmetric faces

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: TRIANGLES_Hmin
REAL (KIND=DBL), SAVE :: Hmin  !minmum height

CONTAINS
	SUBROUTINE Init_Spatial_Grid ()
	IMPLICIT NONE
	CHARACTER (LEN=100) :: String
	INTEGER :: I, J, K, Index, In_TEMP, ios, n1,n2, n_SLAVE, n_MASTER
	INTEGER, DIMENSION(4) :: nd
	REAL (KIND=DBL) :: TMP, xoff, yoff, sl, hl
	REAL (KIND=DBL), DIMENSION(3) :: hlmin
	REAL (KIND=DBL), DIMENSION(4) :: nx,ny
	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: FACES_TMP !temporary arrays: (:,1:2)-node index, (:,3)BCindex, (:,4)periodic pair
	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ELEMENTS ! temporary arrays: dimension (N_ELES, 5): (:,1)-type, (:,2)-PhyID, (:,3:5)-node index
	INTEGER, ALLOCATABLE, DIMENSION (:) :: ID_SLAVE, ID_MASTER ! temporary arrays: faces ID of periodic slave, faces ID of periodic master
	LOGICAL :: INQUE

	OPEN (UNIT=11, FILE=TRIM(FNAME_MSH),STATUS='OLD',ACTION='READ',IOSTAT=ios)

	WRITE (*,*) '**** Initial Spatial Mesh ****'
	DO
		READ (11,*) String
		IF (TRIM(String)=='$Nodes') EXIT
	END DO
	!Read Nodes
	READ (11,*) N_NDS
	ALLOCATE (NODES(N_NDS,2))
	NODES = 0.d0
	DO I = 1, N_NDS
		READ (11,*)Index, NODES(I,1), NODES(I,2), TMP
	END DO

	DO
		READ(11,*)String
		IF (TRIM(String)=='$Elements') EXIT
	END DO
	!Read Elements
	READ (11,*) N_ELES
	ALLOCATE (ELEMENTS(N_ELES,5))
	ALLOCATE (FACES_TMP(N_ELES*3,4)) !(:,1:2)-node index, (:,3)BC flag, (:,4)periodic pair
	FACES_TMP = 0

	ELEMENTS = 0
	N_TRIS = 0
	N_FCS = 0
	DO I = 1, N_ELES
		READ (11,*)Index, ELEMENTS(I,1), In_TEMP, ELEMENTS(I,2)
		BACKSPACE (UNIT=11)
		IF (ELEMENTS(I,1).EQ.1) THEN
			!!only faces at boundaries!!
			READ(11,*) Index, In_TEMP, In_TEMP, In_TEMP, In_TEMP, Elements(I,3:4)
			N_FCS = N_FCS + 1
			FACES_TMP(N_FCS,1:2)=ELEMENTS(I,3:4)
			DO J = 1, NBC
				IF (ELEMENTS(I,2).EQ.BC_PHYID(J)) THEN
					FACES_TMP(N_FCS,3) = J ! BC index
					CYCLE
				END IF
			END DO
		END IF
		IF (ELEMENTS(I,1).EQ.2) THEN
			!triangles
			READ(11,*) Index, In_TEMP, In_TEMP, In_TEMP, In_TEMP, Elements(I,3:5)
			N_TRIS = N_TRIS + 1
			nd(1:3) = ELEMENTS(I,3:5)
			nd(4) = nd(1)
			!extract interior faces
			DO J = 1,3
				INQUE = .FALSE.
				DO K = 1,N_FCS
					n1 = FACES_TMP(K,1)
					n2 = FACES_TMP(K,2)
					IF ((n1.EQ.nd(J)).AND.(n2.EQ.nd(J+1))) THEN
						INQUE = .TRUE.
						CYCLE
					END IF
					IF ((n1.EQ.nd(J+1)).AND.(n2.EQ.nd(J))) THEN
						INQUE = .TRUE.
						CYCLE
					END IF
				END DO
				IF (.NOT.INQUE) THEN
					N_FCS = N_FCS + 1
					FACES_TMP(N_FCS,1) = nd(J)
					FACES_TMP(N_FCS,2) = nd(J+1)
				END IF
			END DO
		END IF
	END DO
	CLOSE(11)

	!FIND periodic pair
	IF (BCID_MASTER > 0) THEN
		n_SLAVE = 0
		n_MASTER = 0
		ALLOCATE (ID_SLAVE(N_FCS),ID_MASTER(N_FCS))

		DO I = 1, N_FCS
			IF (FACES_TMP(I,3).EQ.BCID_MASTER) THEN
				n_MASTER = n_MASTER + 1
				ID_MASTER (n_MASTER) = I
			END IF
			IF (FACES_TMP(I,3).EQ.BCID_SLAVE) THEN
				n_SLAVE = n_SLAVE + 1
				ID_SLAVE (n_SLAVE) = I
			END IF
		END DO

		!pair
		DO I = 1, n_MASTER
			nd(1) = FACES_TMP(ID_MASTER(I),1)
			nd(2) = FACES_TMP(ID_MASTER(I),2)


			xoff = BC_XOFF(BCID_MASTER)
			yoff = BC_YOFF(BCID_MASTER)

			nx(1) = NODES(nd(1),1) + xoff
			ny(1) = NODES(nd(1),2) + yoff
			nx(2) = NODES(nd(2),1) + xoff
			ny(2) = NODES(nd(2),2) + yoff

			DO J = 1, n_SLAVE
				nd(3) = FACES_TMP(ID_SLAVE(J),1)
				nd(4) = FACES_TMP(ID_SLAVE(J),2)
				nx(3) = NODES(nd(3),1)
				ny(3) = NODES(nd(3),2)
				nx(4) = NODES(nd(4),1)
				ny(4) = NODES(nd(4),2)

				IF ((nx(1).EQ.nx(3)).AND.(ny(1).EQ.ny(3))) THEN
					IF ((nx(2).EQ.nx(4)).AND.(ny(2).EQ.ny(4))) THEN
						FACES_TMP(ID_MASTER(I),4)=ID_SLAVE(J)
						FACES_TMP(ID_SLAVE(J),4)=ID_MASTER(I)
						CYCLE
					END IF
				END IF

				IF ((nx(1).EQ.nx(4)).AND.(ny(1).EQ.ny(4))) THEN
					IF ((nx(2).EQ.nx(3)).AND.(ny(2).EQ.ny(3))) THEN
						FACES_TMP(ID_MASTER(I),4)=ID_SLAVE(J)
						FACES_TMP(ID_SLAVE(J),4)=ID_MASTER(I)
						CYCLE
					END IF
				END IF

			END DO
		END DO
		DEALLOCATE (ID_MASTER,ID_SLAVE,STAT=ios)
	END IF

	!Extract triangle information from elements and nodes
	ALLOCATE (TRIANGLES_TAG(N_TRIS,6),TRIANGLES_INF(N_TRIS,7))
	ALLOCATE (TRIANGLES_Hmin(N_TRIS))
	TRIANGLES_TAG = 0
	TRIANGLES_INF = 0.d0
	TRIANGLES_Hmin = 0.d0
	Index = 1
	Hmin = 1.d0
	DO I = 1, N_ELES
		IF (ELEMENTS(I,1).EQ.2) THEN
			!Nodes index
			TRIANGLES_TAG (Index, 1:3) = ELEMENTS (I,3:5)

			nx(1) = NODES(TRIANGLES_TAG(Index,1),1)
			ny(1) = NODES(TRIANGLES_TAG(Index,1),2)
			nx(2) = NODES(TRIANGLES_TAG(Index,2),1)
			ny(2) = NODES(TRIANGLES_TAG(Index,2),2)
			nx(3) = NODES(TRIANGLES_TAG(Index,3),1)
			ny(3) = NODES(TRIANGLES_TAG(Index,3),2)
			nx(4) = nx(1)
			ny(4) = ny(1)
			TRIANGLES_INF(Index,1) = 0.5*(nx(3)*ny(1)-nx(2)*ny(1)+nx(1)*ny(2)-nx(3)*ny(2)-nx(1)*ny(3)+nx(2)*ny(3)) !area
			IF (TRIANGLES_INF(Index,1).lt.0) WRITE (*,*)'Warning: triangle with proper node ordering'

			!outward normal vector
			DO J=1,3
				TRIANGLES_INF(Index,J*2)   = (ny(J+1)-ny(J))/SQRT((nx(J)-nx(J+1))*(nx(J)-nx(J+1))+(ny(J)-ny(J+1))*(ny(J)-ny(J+1)))
				TRIANGLES_INF(Index,J*2+1) = (nx(J)-nx(J+1))/SQRT((nx(J)-nx(J+1))*(nx(J)-nx(J+1))+(ny(J)-ny(J+1))*(ny(J)-ny(J+1)))
			
				sl = SQRT((nx(J)-nx(J+1))**2+(ny(J)-ny(J+1))**2)
				hl = 2.d0*TRIANGLES_INF(Index,1)/sl
				hlmin(J) = hl
				IF (Hmin.GT.hl) Hmin = hl
			END DO
			TRIANGLES_Hmin(Index) = MINVAL(hlmin)

			Index = Index + 1
		END IF
	END DO

	!Extract face information
	ALLOCATE (FACES_TAG(N_FCS,8),FACES_INF(N_FCS))
	FACES_TAG = 0
	FACES_INF = 0.d0
	DO I = 1, N_FCS
		FACES_TAG(I,1:4) = FACES_TMP(I,:)
		nx(1) = NODES(FACES_TAG(I,1),1)
		ny(1) = NODES(FACES_TAG(I,1),2)
		nx(2) = NODES(FACES_TAG(I,2),1)
		ny(2) = NODES(FACES_TAG(I,2),2)
		FACES_INF(I) = SQRT((nx(1)-nx(2))*(nx(1)-nx(2))+(ny(1)-ny(2))*(ny(1)-ny(2)))
	END DO

	!Connect faces and triangle edges
	DO I = 1, N_TRIS
		nd(1) = TRIANGLES_TAG(I,1)
		nd(2) = TRIANGLES_TAG(I,2)
		nd(3) = TRIANGLES_TAG(I,3)
		nd(4) = nd(1)
		DO J = 1, 3
			DO K = 1, N_FCS
				n1 = FACES_TAG(K,1)
				n2 = FACES_TAG(K,2)

				IF ((nd(J).EQ.n1).AND.(nd(J+1).EQ.n2)) THEN
					TRIANGLES_TAG(I,3+J) = K
					FACES_TAG (K,5) = I
					FACES_TAG (K,7) = J
					CYCLE
				END IF

				IF ((nd(J).EQ.n2).AND.(nd(J+1).EQ.n1)) THEN
					TRIANGLES_TAG(I,3+J) = K
					FACES_TAG (K,6) = I
					FACES_TAG (K,8) = J
					CYCLE
				END IF
			END DO
		END DO
	END DO

	!count the number of wall boundary faces and symmetric boundary faces
	N_FCS_WALL_SYM = 0
	DO I = 1, N_FCS
		Index = FACES_TAG(I,3) !BC index
		IF (Index.GT.0) THEN
			IF ((BC_TYP(Index).EQ.1).OR.(BC_TYP(Index).EQ.4)) THEN
				! TYPE = 1: wall; TYPE = 4: symmetry
				N_FCS_WALL_SYM = N_FCS_WALL_SYM + 1
				FACES_TAG(I,4) = N_FCS_WALL_SYM
			END IF
		END IF
	END DO

	N_FCS_I = 0
	N_FCS_B = 0
	DO I = 1, N_FCS
		IF (FACES_TAG(I,3).GT.0) THEN
			N_FCS_B = N_FCS_B + 1
		END IF
	END DO
	N_FCS_I = N_FCS - N_FCS_B


	WRITE (*,100) 'Total 1D faces: ', N_FCS
	100 FORMAT (1X,A20,I7)
	!WRITE (*,101) 'FACE ID', 'Node #1', 'Node #2', 'BC', 'Pair/SYMID', 'TRI(+)', 'TRI(-)', 'LFC(+)', 'LFC(-)', 'LENGTH'
	!101 FORMAT (1X,A8,2A10,A7,A11,4A11,A11)
	!DO I = 1, N_FCS
	!	WRITE (*,102) I, FACES_TAG(I,:), FACES_INF(I)
	!	102 FORMAT (1X,I8,2I10,I7,I11,4I11,ES11.3)
	!END DO
	WRITE (*,*)

	WRITE (*,100) 'Boundary faces: ', N_FCS_B
	WRITE (*,100) 'Interior faces: ', N_FCS_I
	WRITE(*,*)

	WRITE (*,100) 'Total 2D triangles: ', N_TRIS
	!WRITE (*,103) 'TRIANGE ID', 'Node #1', 'Node #2', 'Node #3', 'FACE #1', 'FACE #2', 'FACE #3'
	!103 FORMAT (1X,7A10)
	!DO I = 1, N_TRIS
	!	WRITE(*,104)I,TRIANGLES_TAG(I,:)
	!	104 FORMAT (1X,7I10)
	!END DO
	!WRITE (*,105) 'TRIANGE ID', 'NORMAL VEC (FC#1)', 'NORMAL VEC (FC#2)', 'NORMAL VEC (FC#3)', 'AREA'
	!105 FORMAT (1X,A10,3A25,A11)
	!DO I = 1, N_TRIS
	!	WRITE (*,106)I,'(',TRIANGLES_INF(I,2),',',TRIANGLES_INF(I,3),')','(',TRIANGLES_INF(I,4),',',TRIANGLES_INF(I,5),')' &
	!	,'(',TRIANGLES_INF(I,6),',',TRIANGLES_INF(I,7),')',TRIANGLES_INF(I,1)
	!	106 FORMAT (1X,I10,3(A3,ES10.3,A1,ES10.3,A1),ES11.3)
	!END DO
	WRITE (*,*)
	WRITE(*,108) 'minimun Height: ', Hmin
	108 FORMAT (1X, A16, ES10.2)


	DEALLOCATE(FACES_TMP,STAT=ios)
	DEALLOCATE(ELEMENTS,STAT=ios)

	!DO I = 1, N_TRIS
	!	WRITE(*,*)I, TRIANGLES_Hmin(I)
	!END DO
	WRITE(*,*)MAXVAL(TRIANGLES_Hmin),MINVAL(TRIANGLES_Hmin)
	!STOP
	END SUBROUTINE Init_Spatial_Grid

END MODULE SPATIAL_GRID