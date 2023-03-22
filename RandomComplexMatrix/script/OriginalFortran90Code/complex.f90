!Definition of the MODULE
MODULE complexmatrix
  IMPLICIT NONE
  !definition of the TYPE
  TYPE dmatrix
     INTEGER*4 :: dim(2)
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: mat1
     COMPLEX*16 :: trace
  END TYPE dmatrix
  !FUNCTION
  CONTAINS
  !FUNCTION for RANDOM MATRICES
  TYPE(dmatrix) FUNCTION randomcmatrix(row,col)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: row,col
    INTEGER*4 :: i,j
    REAL*4, DIMENSION(:,:), ALLOCATABLE:: reals, immag

    ALLOCATE(reals(row, col), immag(row,col))
    ALLOCATE(randomcmatrix%mat1(row, col))
    randomcmatrix%dim(1)=row
    randomcmatrix%dim(2)=col
    CALL random_number(reals)
    CALL random_number(immag)
    DO i=1,row
       DO j=1,col
          randomcmatrix%mat1(i,j)=CMPLX(reals(i,j), immag(i,j))
       ENDDO
    ENDDO
    DEALLOCATE(reals, immag)
  END FUNCTION randomcmatrix
  !FUNCTION for HAND-WRITTEN MATRICES
  TYPE(dmatrix)  FUNCTION handcmatrix(row,col)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: row,col
    INTEGER :: i,j
    REAL*4 :: reals, immag
    
    ALLOCATE(handcmatrix%mat1(row,col))
    handcmatrix%dim(1)=row
    handcmatrix%dim(2)=col
    DO i=1,row
       DO j=1,col
          PRINT*, 'Insert the REAL part of the element (', i, ',',j,'):'
          READ*, reals
          PRINT*, 'Insert the IMAGINARY part of the element (', i, ',',j,'):'
          READ*, immag
          handcmatrix%mat1(i,j) = CMPLX(reals, immag)
       ENDDO
    ENDDO 
  END FUNCTION handcmatrix
  !SUBROUTINE for TRACE
  SUBROUTINE tracecmatrix(matr2)
    IMPLICIT NONE
    TYPE(dmatrix) :: matr2
    INTEGER*4 :: i
    COMPLEX*16 :: tr

    tr = CMPLX(0.0,0.0)
    DO i=1, matr2%dim(1)
        tr = tr + matr2%mat1(i,i)
    ENDDO
    matr2%trace = tr  
  END SUBROUTINE tracecmatrix
  !FUNCTION for the ADJOINT
  TYPE(dmatrix) FUNCTION adjointcmatrix (mat2)
    IMPLICIT NONE
    INTEGER*4 :: i,j
    TYPE(dmatrix), INTENT(IN) :: mat2
    
    adjointcmatrix%dim(1) = mat2%dim(2)
    adjointcmatrix%dim(2) = mat2%dim(1)
    ALLOCATE(adjointcmatrix%mat1(adjointcmatrix%dim(1), adjointcmatrix%dim(2)))
    DO i=1, adjointcmatrix%dim(1)
       DO j=1, adjointcmatrix%dim(2)
          adjointcmatrix%mat1(i,j) = CMPLX(REAL(mat2%mat1(j,i)), -AIMAG(mat2%mat1(j,i)))
       ENDDO
    ENDDO
  END FUNCTION adjointcmatrix
END MODULE complexmatrix


!MAIN PROGRAM
PROGRAM interface
  USE complexmatrix
  IMPLICIT NONE
  !INTERFACE RANDOM MATRIX
  INTERFACE OPERATOR (.randomcmatrix.)
     MODULE PROCEDURE randomcmatrix
  END INTERFACE OPERATOR (.randomcmatrix.)
  !INTERFACE HAND-WRITTEN MATRIX
  INTERFACE OPERATOR (.handcmatrix.)
     MODULE PROCEDURE handcmatrix
  END INTERFACE OPERATOR (.handcmatrix.)
  !INTERFACE TRACE
  INTERFACE OPERATOR 
     MODULE PROCEDURE tracecmatrix
  END INTERFACE OPERATOR
  !INTERFACE ADJOINT
  INTERFACE OPERATOR (.adjointcmatrix.)
     MODULE PROCEDURE adjointcmatrix
  END INTERFACE OPERATOR (.adjointcmatrix.)
  
  INTEGER*4 :: row, col
  TYPE(dmatrix) :: matrice1, adjmatrice1
  INTEGER :: cho, i ,j
  CHARACTER*40 :: test='results.txt'

  PRINT*, 'Please, inster the number of rows for the matrix: '
  READ*, row
  PRINT*, 'Please, inster the number of columns for the matrix: '
  READ*, col
  PRINT*, 'What do you want to do?'
  PRINT*, '1) RANDOM COMPLEX MATRIX'
  PRINT*, '2) HAND-WRITTEN COMPLEX MATRIX'
  READ*, cho
  IF (cho==1) THEN
     matrice1 = row.randomcmatrix.col
  ELSEIF (cho==2) THEN
     matrice1 = row.handcmatrix.col
  ELSE
     PRINT*, 'You make no choice'
  ENDIF
  CALL tracecmatrix(matrice1)
  PRINT*, 'Your Matrix is'
  DO i=LBOUND(matrice1%mat1,1), UBOUND(matrice1%mat1,1)
        WRITE(*,*) (matrice1%mat1(i,j), j=LBOUND(matrice1%mat1,2), UBOUND(matrice1%mat1,2))
  ENDDO
  PRINT*, 'The trace of this matrix is: ', matrice1%trace
  adjmatrice1 = .adjointcmatrix.matrice1
  CALL tracecmatrix(adjmatrice1)
  PRINT*, 'Your ADJOINT Matrix is'
  DO i=LBOUND(adjmatrice1%mat1,1), UBOUND(adjmatrice1%mat1,1)
       WRITE(*,*) (adjmatrice1%mat1(i,j), j=LBOUND(adjmatrice1%mat1,2), UBOUND(adjmatrice1%mat1,2))
  ENDDO
  PRINT*, 'The trace of this adjoint matrix is: ', adjmatrice1%trace
  PRINT*, ''
  PRINT*,'Do you wanna print it in a TXT file? (1) for yes'
  READ*, cho
  IF (cho==1) THEN
     test = TRIM(test)
     OPEN(UNIT=50, file=test, status='unknown')
     WRITE(50,*) 'The matrix is:'
     DO i=LBOUND(matrice1%mat1,1), UBOUND(matrice1%mat1,1)
        WRITE(50,*) (matrice1%mat1(i,j), j=LBOUND(matrice1%mat1,2), UBOUND(matrice1%mat1,2))
     ENDDO
     WRITE(50,*) ''
     WRITE(50,*) 'The dimensions of the matrix is:', matrice1%dim(1),'x', matrice1%dim(2)
     WRITE(50,*) 'The trace of the matrix is: ', matrice1%trace
     WRITE(50,*) ''
     WRITE(50,*) 'The ADJOINT Matrix is:'
     DO i=LBOUND(adjmatrice1%mat1,1), UBOUND(adjmatrice1%mat1,1)
        WRITE(50,*) (adjmatrice1%mat1(i,j), j=LBOUND(adjmatrice1%mat1,2), UBOUND(adjmatrice1%mat1,2))
     ENDDO     
     WRITE(50,*) ''
     WRITE(50,*) 'The dimensions of the matrix is:', adjmatrice1%dim(1),'x', adjmatrice1%dim(2)
     WRITE(50,*) 'The trace of the matrix is: ', adjmatrice1%trace
     CLOSE(50)
  ENDIF
END PROGRAM interface


