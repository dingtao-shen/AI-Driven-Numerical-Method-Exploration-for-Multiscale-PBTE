MODULE NODAL_FUNCTION
!Define the basis function
USE CONSTANT
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
IMPLICIT NONE

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: NODFUN_FC_REF ! coefficients of face nodal function in reference line
                                                                    ! -1 -------- 1
                                                                    ! DIMENSION(NDOF_FC,NDOF_FC)
                                                                    ! NODFUN_FC_REF(:,1) + NODFUN_FC_REF(:,2)*t + ... +NODFUN_FC_REF(:,n)*t**(n-1)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: NODE_FC_REF !nodal points in reference line
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:),SAVE :: NODFUN_TRI_REF !coefficients of element nodal function in reference triangle
                                                                    ! (0,1)
                                                                    !   |\
                                                                    !   | \
                                                                    !   |  \
                                                                    !   |___\ (1,0)
                                                                    !  (0,0)
                                                                    ! DIMESION (NDOF_TRI,NDOF_TRI)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: NODE_TRI_REF !nodal points in reference triangle, DIMESION(NDOF_TRI,2)

CONTAINS

SUBROUTINE Init_Basis_Function ()
IMPLICIT NONE
!INTEGER :: I, J, K, L, M

WRITE (*,*)'**** Initial Basis Function ****'
ALLOCATE (NODFUN_FC_REF(NDOF_FC,NDOF_FC),NODE_FC_REF(NDOF_FC))
ALLOCATE (NODFUN_TRI_REF(NDOF_TRI,NDOF_TRI),NODE_TRI_REF(NDOF_TRI,2))

NODFUN_FC_REF = 0.d0
NODE_FC_REF = 0.d0
NODFUN_TRI_REF = 0.d0
NODE_TRI_REF = 0.d0

!Nodal function in reference line
CALL Basis_Function_FC (NDOF_FC,NODE_FC_REF,NODFUN_FC_REF)
!Nodal function in reference triangle
CALL Basis_Function_TRI(DEG,NDOF_TRI,NODE_TRI_REF,NODFUN_TRI_REF)


END SUBROUTINE Init_Basis_Function

SUBROUTINE Basis_Function_FC (NN,NODE_REF,NODFUN_REF)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NN
REAL (KIND=DBL), INTENT(OUT), DIMENSION(NN,NN) :: NODFUN_REF
REAL (KIND=DBL), INTENT(OUT), DIMENSION(NN) :: NODE_REF
INTEGER :: I,J,K,L
REAL (KIND=DBL) :: b0, b1
REAL (KIND=DBL), DIMENSION(NN) :: a

NODE_REF = 0.d0
NODFUN_REF = 0.d0

!nodes
DO I = 1, NN
	NODE_REF(I) = -1.d0 + 2.d0*REAL(I-1,DBL)/REAL(NN-1,DBL)
END DO

!function
DO I = 1, NN
	NODFUN_REF(I,1) = 1.d0
	K = 0
	DO J = 1, NN
		IF (J.NE.I) THEN
			a(1:K+1) = NODFUN_REF(I,1:K+1)
			b0 = REAL(2*J-NN-1,DBL)/REAL(2*(J-I),DBL)
			b1 = -REAL(NN-1,DBL)/REAL(2*(J-I),DBL)

			NODFUN_REF(I,1) = a(1)*b0
			NODFUN_REF(I,K+2) = a(K+1)*b1

			DO L = 2,K+1
				NODFUN_REF(I,L) = a(L)*b0 + a(L-1)*b1
			END DO
			K = K + 1
		END IF
	END DO
END DO

END SUBROUTINE

SUBROUTINE Basis_Function_TRI (RR,NN,NODE_REF,NODFUN_REF)
IMPLICIT NONE
INTEGER, INTENT(IN) :: RR,NN
REAL (KIND=DBL), INTENT(OUT), DIMENSION(NN,NN) :: NODFUN_REF
REAL (KIND=DBL), INTENT(OUT), DIMENSION(NN,2) :: NODE_REF
INTEGER :: I,J,K,count,L,M,N,pos,C,pos1,pos2
REAL (KIND=DBL), DIMENSION(NN) :: a, b, tmp
REAL (KIND=DBL) :: d0,d1,d2

NODE_REF = 0.d0
NODFUN_REF = 0.d0

count = 1
DO J = 0, RR
	DO I = 0, RR-J
		NODE_REF(count,1) = REAL(I,DBL)/REAL(RR,DBL)
		NODE_REF(count,2) = REAL(J,DBL)/REAL(RR,DBL)
		count = count + 1
	END DO
END DO

count = 1
DO J = 0, RR
	DO I = 0, RR-J
		K = RR - I - J
		a = 0.d0
		b = 0.d0
		c = 0.d0

		IF (I.NE.0) THEN
			a(2) = REAL(RR,DBL)/REAL(I,DBL)
			IF (I.GT.1) THEN
				DO L = 1, I-1
					tmp = a
					d0 = REAL(L,DBL)/REAL(L-I,DBL)
					d1 = -REAL(RR,DBL)/REAL(L-I,DBL)

					a(1) = tmp(1)*d0
					a((L+1)*(L+2)/2+1) = tmp(L*(L+1)/2+1)*d1

					DO M = 1, L
						a(M*(M+1)/2+1) = tmp(M*(M+1)/2+1)*d0 + tmp((M-1)*M/2+1)*d1 
					END DO
				END DO
			END IF
		ELSE
			a(1) = 1.d0
		END IF

		IF (J.NE.0) THEN
			b(3) = REAL(RR,DBL)/REAL(J,DBL)
			IF (J.GT.1) THEN
				DO M = 1, J-1
					tmp = b
					d0 = REAL(M,DBL)/REAL(M-J,DBL)
					d1 = -REAL(RR,DBL)/REAL(M-J,DBL)

					b(1) = tmp(1)*d0
					b((M+2)*(M+3)/2) = tmp((M+1)*(M+2)/2)*d1

					DO L = 1, M
						b((L+1)*(L+2)/2) = tmp((L+1)*(L+2)/2)*d0 + tmp(L*(L+1)/2)*d1
					END DO
				END DO
			END IF
		ELSE
			b(1) = 1.d0
		END IF

		DO L = 0, I
			DO M = 0, J
				pos = (L+M)*(L+M+1)/2+M+1
				NODFUN_REF(count,pos) = NODFUN_REF(count,pos) + &
					a(L*(L+1)/2+1)*b((M+1)*(M+2)/2)
			END DO
		END DO


		IF (K.NE.0) THEN
			DO N = 0, K-1
				tmp = NODFUN_REF(count,:)
				d0 = REAL(N-RR,DBL)/REAL(N-K,DBL)
				d1 = REAL(RR,DBL)/REAL(N-K,DBL)
				d2 = REAL(RR,DBL)/REAL(N-K,DBL)

				DO C = 0, RR
					DO M = 0, C
						L = C - M
						pos = (L+M)*(L+M+1)/2+M+1
						pos1 = (L-1+M)*(L-1+M+1)/2+M+1
						pos2 = (L+M-1)*(L+M)/2+M

						NODFUN_REF(count,pos)= tmp(pos)*d0
						IF (L.GT.0) NODFUN_REF(count,pos) = NODFUN_REF(count,pos) + tmp(pos1)*d1
						IF (M.GT.0) NODFUN_REF(count,pos) = NODFUN_REF(count,pos) + tmp(pos2)*d2
					END DO
				END DO

			END DO
		END IF
		count = count + 1
	END DO
END DO

END SUBROUTINE

END MODULE NODAL_FUNCTION