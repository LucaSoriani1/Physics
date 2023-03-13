!======================================================================!
!This program aims to compute the algorithm on the renormalization
!group, namely ''real space renormalization group''. The program use
!this algorithm on a tranverse field Ising Model as exampel of quantum
!many-body problem. The program consider a fixed number of particles,
!then compute the Ising model. After that, the program computes
!an reasonable approximation for the solution of the many body problem.
!======================================================================!



!======================================================================!
!In this module there are the main derived type used in the calculation
!(the ''dmatrix'' type that contains the double complex matrix) and the
!usefull tools for the computation, as the subroutine that calculates
!the trace or the subroutine used to substituite the matrix to each
!other.
!======================================================================!
MODULE tools

  !--------------------------------------------------------------------!
  !Creation of a more elaborated type of variable, as the first
  !exercise asked
  TYPE dmatrix

     !Two integer, respectivelly, for the row and the columns of
     !the given matrix
     INTEGER*8 :: row, col

     !the allocatable double complex matrix
     DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: mat

     !the double precision array to store the eigenvectors
     DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE :: eigval

     !Double complex variable for the determinant
     COMPLEX*16 :: det

     !Double complex variable for the trace
     COMPLEX*16 :: trace

  END TYPE dmatrix
  !--------------------------------------------------------------------!


CONTAINS  !Subroutines used in the exercises

  !=====================================================================!
  !This subroutine checks if a condition is true or false. As input it
  !take the condition and a message, that is the error message. In the
  !main program, the condition is represented by the wrong condition: if
  !this wrong condition is true, the program will print the error message
  !given as input and will print the position of the check in order to
  !simplify the search of the problem.
  !=====================================================================!
  SUBROUTINE checkpoint(condition, error)

    IMPLICIT NONE

    !This boolean variable is defined by input in the main program.
    LOGICAL, INTENT(IN) :: condition

    !This integer is necessary to the refer the position of the checkpoint.
    !the SAVE permits to save of the variable after the call in the main
    !program.
    INTEGER, SAVE :: pos=1

    !Set of character for the error message, represented by a 
    !character of unspecified length.
    CHARACTER(*), INTENT(IN) :: error

    !Checkpoint statement
    IF(condition .EQV. .TRUE.) THEN
       PRINT*,'ERROR: ',  error
       WRITE(*, '("At position: ",I0)') pos
       STOP
    ENDIF

    !Increase the position of the error
    pos = pos + 1
  END SUBROUTINE checkpoint
  !=====================================================================!

  !=====================================================================!
  !This subroutine is just to keep clean the main program. This creates
  !the main matrices: the Pauli's matrices and the identity matrix.
  !It takes as input and output four 'dmatrix' (type created before)
  !and an integer that will be the dimension of the matrices.
  !=====================================================================!
  SUBROUTINE pauli(x, z, id, dim)

    !Four 'dmatrix' input and output
    TYPE(dmatrix), INTENT(INOUT) :: x, z, id

    !Dimension of the matrices
    INTEGER*8 :: dim

    !Integer used for the checks
    INTEGER*8 :: info

    !Initialization of the dimension of the matrices:
    !Sigma-X
    x%row=dim
    x%col=dim
    !Sigma-Z
    z%row=dim
    z%col=dim
    !Identity
    id%row=dim
    id%col=dim

    !ALLOCATION of the matrices
    ALLOCATE(x%mat(dim, dim), z%mat(dim, dim), id%mat(dim,dim), STAT=info)

    !CHECKPOINT: Checking the allocation
    CALL checkpoint(info/=0, "Somenthing wrong with the allocation of the Pauli's &
         & matrices/identity matrix. Maybe not enough memory. Please, check it")

    !Initialization of the matrices

    !Sigma-X
    x%mat(1,1) = (0d0, 0d0)
    x%mat(1,2) = (1d0, 0d0)
    x%mat(2,1) = (1d0, 0d0)
    x%mat(2,2) = (0d0, 0d0)

    !Sigma-Z
    z%mat(1,1) = (1d0, 0d0)
    z%mat(2,1) = (0d0, 0d0)
    z%mat(1,2) = (0d0, 0d0)
    z%mat(2,2) = (-1d0, 0d0)

    !Identity
    id%mat(1,1) = (1d0, 0d0)
    id%mat(2,1) = (0d0, 0d0)
    id%mat(1,2) = (0d0, 0d0)
    id%mat(2,2) = (1d0, 0d0)

    !Returning the result
    RETURN

  END SUBROUTINE pauli
  !=====================================================================!


  !=====================================================================!
  !This subroutine is made to simplify the creation of the Z-part and
  !the X-part of the Hamiltonian operator. It replaces the first input
  !with the second. Because the first input could have a lower dimension
  !with respect the second, it re-creates the matrix and copy the second
  !input in the first. As input it takes two 'dmatrix' and gives back
  !the first input equal to the second and the second is deallocated.
  !=====================================================================!
  SUBROUTINE substi(mat1, mat2)

    IMPLICIT NONE

    !Two 'dmatrx' as input and output
    TYPE(dmatrix), INTENT(INOUT) :: mat1, mat2

    !Integer used to check the allocation
    INTEGER*4 :: info

    !Deallocation of the first input matrix
    DEALLOCATE(mat1%mat, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Somenthing wrong with the deallocation in ''substi''. Please check it')

    !Defining the dimensions of the new matrix as the dimensions of
    !the second matrix
    mat1%row=mat2%row
    mat1%col=mat2%col

    !Allocation of the new matrix with the new dimension
    ALLOCATE(mat1%mat(mat1%row, mat1%col), STAT=info)

    !CHECKPOINT: checking the allocation of the copied matrix
    CALL checkpoint(info/=0, 'Somenthing wrong with the allocation in ''substi''. Please check it')

    !Copying the second matrix in the first
    mat1%mat = mat2%mat

    !Deallocation of the second matrix
    DEALLOCATE(mat2%mat, STAT=info)

    !CHECKPOINT: checking the deallocation of the copied matrix
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''substi''. Please, check it')

    !Returning the result
    RETURN

  END SUBROUTINE substi
  !=====================================================================!


  !=====================================================================!
  !This subroutine computes the tensor product between two matrices.
  !It take as input two 'dmatrix' types and gives back a 'dmatrix'
  !that contains the matrix as tensor product of the first two.
  !=====================================================================!
  SUBROUTINE t_prod(mat1, id, tempor)

    IMPLICIT NONE

    !Inputs that I want multiply tensorially
    TYPE(dmatrix), INTENT(IN) :: mat1, id

    !Output of the tensor product
    TYPE(dmatrix), INTENT(OUT) :: tempor

    !Integer used for the DO cycles
    INTEGER*8 :: ii, jj, kk, ll

    !Integer used to check the allocation
    INTEGER*4 :: info

    !Defining the dimension of the output matrix as the product of the
    !dimension of the two input matrices
    tempor%row=INT(mat1%row*id%row)
    tempor%col=INT(mat1%col*id%col)

    !Allocation of the outmatrix
    ALLOCATE(tempor%mat(tempor%row, tempor%col), STAT=info)

    tempor%mat=CMPLX(0d0, 0d0)

    !CHECKPOINT: checking the allocation of the matrix
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''t_prod''. Please check it1')
    !Tensor product
    DO ii=1, mat1%row
       DO jj=1, mat1%col
          DO kk=1, id%row
             DO ll=1, id%col
                tempor%mat(kk+(ii-1)*(id%row), ll+(jj-1)*(id%col)) = mat1%mat(ii,jj)*id%mat(kk, ll)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !Returning the results
    RETURN

  END SUBROUTINE t_prod
  !=====================================================================!


  !=====================================================================!
  !This subroutine compute the trace of a input matrix. It takes as
  !input the corrisponding 'dmatrix' type and compute the trace and
  !store it inside the type.
  !=====================================================================!
  SUBROUTINE trace(matrice)

    !Input and output: the matrix which the program computes the trace
    TYPE(dmatrix), INTENT(INOUT) :: matrice

    !Double precision real to store the trace
    REAL*8 :: traccia

    !Integer used in the DO cycels
    INTEGER*8 :: ii

    !Initialization of the trace to zero
    traccia = 0d0

    !Computation of the trace
    DO ii=1, matrice%row
       traccia = traccia +matrice%mat(ii, ii)
    ENDDO

    !Storing the trace in the type
    matrice%trace = traccia

    !Returning the result
    RETURN

  END SUBROUTINE trace
  !=====================================================================!


  !====================================================================!
  !This subroutine copy the input matrix in another one. The suborutine
  !takes an input a 'dmatrix' type. As output gives the same 'dmatrix'
  !but copied in another one.
  !====================================================================!
  SUBROUTINE copy(matA, matB)

    IMPLICIT NONE

    !Input: 'dmatrix' to copy
    TYPE(dmatrix), INTENT(IN) :: matA

    !Input and Output: copied 'dmatrix'
    TYPE(dmatrix), INTENT(OUT) :: matB

    !Integer used for the checkpoint
    INTEGER*8 :: info

    !Copying the dimension
    matB%row=matA%row
    matB%col=matA%col

    !Allocation of the copied matrix
    ALLOCATE(matB%mat(matB%row, matB%col), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''copy''. Please, check it.')

    !Copying the matrix
    matB%mat = matA%mat

    !Returing the result
    RETURN

  END SUBROUTINE copy
  !====================================================================!


  !====================================================================!
  !This subroutine takes a matrix making it the hermitian transpose. It
  !takes as input the input and output the matrix, that will
  !transformed
  !====================================================================!
  SUBROUTINE hermi(mat1)

    IMPLICIT NONE

    !Input and output: the matrix to make hermitian
    TYPE(dmatrix), INTENT(INOUT) :: mat1

    !Temporary matrix
    TYPE(dmatrix) :: tempo

    !Integer for the DO cycle
    INTEGER*8 :: ii, jj, info

    !Integer for the checkpoint

    !Intering the rows and columns
    tempo%row=mat1%col
    tempo%col=mat1%row

    !Allocation of the temporary matrix
    ALLOCATE(tempo%mat(tempo%row, tempo%col), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''hermi''. Please, check it.')

    !Making the temporary matrix hermitian
    DO ii=1, tempo%row
       DO jj=1, tempo%col
          tempo%mat(ii, jj) = CONJG(mat1%mat(jj, ii))
       ENDDO
    ENDDO

    !Deallocation of the input
    DEALLOCATE(mat1%mat, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''hermi''. Please, check it.')

    !Defining the input as the temporary matrix
    CALL copy(tempo, mat1)

    !Deallocation of the temporary matrix
    DEALLOCATE(tempo%mat, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''hermi''. Please, check it.')

    !Returning the result
    RETURN

  END SUBROUTINE hermi
  !====================================================================!


  !====================================================================!
  !This subroutine solves the eigenproblem. It takes as input the
  !'dmatrix' that the user wants to diagonalize, and then the subroutine
  !solves the eigenproblem using the ZHEEV subroutine from LAPACK.
  !====================================================================!
  SUBROUTINE eigenp(mat1)

    IMPLICIT NONE

    !Input and output: the 'dmatrix' type to diagonalize
    TYPE(dmatrix), INTENT(INOUT) :: mat1

    !Variable used for the ZHEEV subroutine
    INTEGER*8 :: lwork=-1, lwmax=2000

    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: work

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rwork

    INTEGER*4 :: info

    !Temporary matrix
    TYPE(dmatrix) :: trix1

    !COpying the input into the temporary matrix
    CALL copy(mat1, trix1)

    !Deallocation of the input matrix
    DEALLOCATE(mat1%mat, STAT=info)

    !CHECKPOINT:checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''eigenp''. Please check it.')
    
    !Defining lwrok=-1 in order to compute the optimal 'lwork'
    !for the computation of the eigenvalues
    lwork=-1

    !Allocation of the array for the computation of the eigenvalues
    ALLOCATE(trix1%eigval(trix1%row), work(lwmax), rwork(3*trix1%row-2), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in the ''eigenp''. Please, check it.')

    !Computing the optimal 'lwork'
    CALL ZHEEV('V', 'U', trix1%col, trix1%mat, trix1%row, trix1%eigval, work, &
         & lwork, rwork, info)

    !CHECKPOINT: checking the work of the ZHEEV subroutine
    CALL checkpoint(info/=0, 'Something wrong with the first ZHEEV. Please, check it')

    !Defining of the optimal 'lwork'
    lwork = MIN(lwmax, INT(work(1)))

    !Deallocation of the 'work' array
    DEALLOCATE(work, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''eigenp''. Please check it.')

    !Allocation of the 'work' array of the right size
    ALLOCATE(work(MAX(1, lwork)), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''eigenp''. Please, check it.')

    !Computing of the eigenvalues
    CALL ZHEEV('V', 'U', trix1%col, trix1%mat, trix1%row, trix1%eigval, work, &
         & lwork, rwork, info)

    !CHECKPOINT: checking the work of the ZHEEV subroutine
    CALL checkpoint(info/=0, 'Something wrong with the second ZHEEV. Please, check it.')

    !Deallocation of the arrays used in the computation
    DEALLOCATE(work, rwork, STAT=info)

    !CHECKPOINT: checking the deallocation of the arrays
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''eigenp''. Please, check it.')

    !copying the result into the output
    CALL copy(trix1, mat1)

    !Allocation of the output eigenvalues Array)
    ALLOCATE(mat1%eigval(mat1%row), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''eigenp''. Please, check it.')

    !Copying the eigenvalues
    mat1%eigval = trix1%eigval

    !Deallocation of the temporary matrix
    DEALLOCATE(trix1%mat, trix1%eigval, STAT=info)
    
    !CHECKPOINT: checking the deallocation of the arrays
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''eigenp''. Please, check it.')
    
    !Returning the results
    RETURN

  END SUBROUTINE eigenp
  !====================================================================!



  !====================================================================!
  !This subroutine executes the matrix product between two matrix. It
  !takes as input the two matix and output returns the product between
  !them. As function, this subroutine will used the MATMUL.
  !====================================================================!
  SUBROUTINE multi(trix1, trix2, resu)

    IMPLICIT NONE

    !INPUT: the two matrices that will be multiply
    TYPE(dmatrix), INTENT(IN) :: trix1, trix2

    !OUTPUT: the result of the multiplication
    TYPE(dmatrix), INTENT(OUT) :: resu

    !Integer used for the checkpoint
    INTEGER*8 :: info

    !CHECKPOINT: checking the right dimension for the MATMUL operation
    CALL checkpoint(SIZE(trix1%mat,2) /= SIZE(trix2%mat,1), 'MATMUL cannot multiply the &
         &matrix beacuse the number of columns of the first matrix is not equal to the&
         & number of rows of the second matrix. Please, check it.')

    !Defining the output
    resu%row=trix1%row
    resu%col=trix2%col

    !Allocation of the output matrix
    ALLOCATE(resu%mat(resu%row, resu%col), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''multi''.&
         & Please, check it.')

    !Multiplication
    resu%mat = MATMUL(trix1%mat, trix2%mat)

    !Returning the result
    RETURN

  END SUBROUTINE multi
  !====================================================================!


  !====================================================================!
  !This subroutine computes the reduced density matrix for a given
  !matrix. As input it takes the matrix that will be reduced, the index
  !of the sub-system, and the dimension of the sub-system. As ouput it
  !returns the reduced density matrix for the chosen sub-system
  !====================================================================!
  SUBROUTINE reduced(rho, redu, indk, dim)

    IMPLICIT NONE

    !Input: the matrix that will be reduced
    TYPE(dmatrix), INTENT(IN) :: rho

    !Output: the reduced matrix
    TYPE(dmatrix), INTENT(OUT) :: redu

    !Input: integer to mark the sub-system
    INTEGER*8, INTENT(IN) :: indk

    !Input: dimension of the sub-system
    INTEGER*8, INTENT(IN) :: dim

    !Integer used for the checkpoint and for the DO cycles
    INTEGER*8 :: info, ii, jj, kk

    !Integer used to identify the rows and the columns
    INTEGER*8 :: row, col

    !Defining the dimension of the reduced density matrix
    redu%row = INT(SQRT(DBLE(SIZE(rho%mat, 1))))
    redu%col = INT(SQRT(DBLE(SIZE(rho%mat, 2))))

    !Allocation of the reduced density matrix
    ALLOCATE(redu%mat(redu%row, redu%col), STAT=info)

    !CHECKPOINT: checking the allocation
    CALL checkpoint(info/=0, 'Something wrong with the allocation in ''reduced''. &
         &Please, check it.')

    !Initialization of the reduced matrix as null matrix
    redu%mat = CMPLX(0d0, 0d0)

    !Computation of the reduced density matrix
    DO ii = 1, SIZE(redu%mat, 1)
       DO jj= 1, SIZE(redu%mat, 2)
          DO kk=1, dim
             !row
             row = 1 + (MOD((ii-1), (dim)**(indk-1))) + (kk-1)*(dim*(indk-1)) + &
                  & dim*(ii-1-(MOD((ii-1), (dim)**(indk-1)))) + (kk-1)*(2-indk)

             !col
             col = 1 + (MOD((jj-1), (dim)**(indk-1))) + (kk-1)*(dim*(indk-1)) + &
                  & dim*(jj-1-(MOD((jj-1), (dim)**(indk-1)))) + (kk-1)*(2-indk)

             !reduced
             redu%mat(ii,jj) = redu%mat(ii,jj) + rho%mat(row, col)
          ENDDO
       ENDDO
    ENDDO

    !Checking the trace of the reduced density matrix
    CALL trace(redu)

    !CHECKPOINT: checking the trace of the reduced density matrix is 1
    CALL checkpoint(ABS(ABS(redu%trace)-1d0) > 0.0004, 'The trace of the reduced &
         &density matrix in ''reduced'' is not 1. Please, check it.')

    !Checking if the reduced density matrix is hermitian
    DO ii=1, redu%row
       DO jj=1, redu%col
          CALL checkpoint(ABS(redu%mat(ii, jj) - CONJG(redu%mat(jj,ii))) > 0.1, &
               &'The reduced density matrix is not hermitian. Please, check it.')
       ENDDO
    ENDDO

    !Returning the result
    RETURN
    
  END SUBROUTINE reduced
  !====================================================================!
  

END MODULE tools
!======================================================================!


!======================================================================!
!This module contains the two subroutines used to compute the Z-part
!and the X-part of the Ising Model.
!======================================================================!
MODULE ising

  !Using the previous module
  USE tools

CONTAINS

  !=====================================================================!
  !This subroutines creates the Z-part of the Ising Model. It takes as
  !input two 'dmatrix', i.e the Sigma-Z Pauli's matrix and the identity
  !matrix used to define the Z-part matrix, the integer ''N'', i.e. the
  !number of particles (thus the iteraction). As output, the subroutine
  !gives back the Z-part matrix.
  !=====================================================================!
  SUBROUTINE make_ham_z(mat1, id, N, out)

    IMPLICIT NONE

    !Input: the Sigma-Z Pauli's matrix and the identity matrix
    TYPE(dmatrix), INTENT(IN) :: mat1, id

    !Input: Number of particles
    INTEGER*8, INTENT(IN) :: N

    !Output: the Z-part of the Ising Model
    TYPE(dmatrix), INTENT(OUT) :: out

    !Integer used to store the dimension of the output matrix
    INTEGER*8 :: dim1, dim2

    !Integer used in the DO cycles
    INTEGER*8 :: ii, jj, kk

    !Temporary matrix used for the computation
    TYPE(dmatrix) :: trix1, trix2

    !Integer used to check the deallocation/allocation
    INTEGER*4 :: info

    !Defining the dimensions of the output matrix
    dim1=INT(mat1%row**N)
    dim2=INT(mat1%col**N)
    out%row=dim1
    out%col=dim2

    !Allocation of the Z-part matrix
    ALLOCATE(out%mat(dim1, dim2), STAT=info)

    !CHECKPOINT: checking the allocation of the matrix
    CALL checkpoint(info/=0, 'Somenthing wrong with the allocation of the &
         & Z-part matrix. Maybe not enough memory. Please, check it.')

    !Initialization of the matrix as null matrix
    out%mat = (0d0, 0d0)

    !Starting the cycle about the computation of the Z-part.
    !I treat three cases: first when the index is ii=1, i.e. I am
    !considering the first particle, the second when the index is
    !ii=2, i.e. I am considering the second particle, and then
    !when the index is ii>2. At the end of each cycle, The program
    !sums the result to the previous one. 
    DO ii=1, N
       !CASE 1: First particle. The first tensor product is between
       !the Sigma-Z Pauli's matrix and the identity. Thus, the program
       !computes the order tensor products between the results and the
       !identity matrix
       IF (ii==1) THEN

          !Tensor product between the Sigma-Z and the identity
          CALL t_prod(mat1, id, trix1)

          IF(N/=2)THEN
             !Other tensor products
             DO jj=1, N-2

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)

             ENDDO
          ENDIF

          !CASE 2: Second particle. The first tensor product is between
          !the identity and Sigma-Z Pauli's matrix. Thus, the program
          !computes the other tensor products between the result and
          !the identity matrix
       ELSEIF (ii==2) THEN

          !Tensor product between the identity and the Sigma-Z
          CALL t_prod(id, mat1, trix1)

          IF(N/=2)THEN
             !Other tensor products
             DO jj=1, N-2

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)
             ENDDO
          ENDIF


          !CASE 3: other particles. The program computer the first ii-2
          !tensor products and then computes the tensor product between
          !the result and the Sigma-Z Pauli's matrix. Then, if I am not
          !considering the last particle, the program computes the tensor
          !products between the result and the other identity matrix.
          !If I am considering the last particle, the program ends
       ELSE

          !Initialize the temporary matrix as the 2X2 identity matrix.
          trix1%row=id%row
          trix1%col=id%col
          ALLOCATE(trix1%mat(trix1%row, trix1%col), STAT=info)

          !CHECKPOINT: checking the allocation
          CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it')

          trix1%mat = id%mat

          !First ii-2 tensor products between the identities
          DO jj=1, ii-2

             CALL t_prod(trix1, id, trix2)

             CALL substi(trix1, trix2)

          ENDDO

          !The tensor product between the first ii-2 identity matrix
          !and the Sigma-Z
          CALL t_prod(trix1, mat1, trix2)

          CALL substi(trix1, trix2)

          !If I am not considering the last particle, the program
          !computes the remaining tensor products
          IF (ii /= N) THEN

             !Other tensor products
             DO jj=1, N-ii

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)
             ENDDO
          ENDIF

       ENDIF

       !Storing the output matrix, aka the Z-part of the Ising Model
       DO jj=1, out%row
          DO kk=1, out%col
             out%mat(jj,kk) = out%mat(jj, kk) + trix1%mat(jj, kk)
          ENDDO
       ENDDO

       !DEALLOCATION of the temporary matrix in order to
       !restart the cycle correctly
       DEALLOCATE(trix1%mat, STAT=info)

       !CHECKPOINT: checking the deallocation
       CALL checkpoint(info/=0, 'Something wrong with the deallocation. &
            &Please, check it.')
    ENDDO

    !Computation of the trace of the output matrix
    CALL trace(out)

    !Returning the result
    RETURN

  END SUBROUTINE make_ham_z
  !=====================================================================!



  !=====================================================================!
  !This subroutines creates the X-part of the Ising Model. It takes as
  !input two 'dmatrix', i.e the Sigma-X Pauli's matrix and the identity
  !matrix used to define the X-part matrix, the integer ''N'', i.e. the
  !number of particles (thus the iteraction). As output, the subroutine
  !gives back the X-part matrix.
  !=====================================================================!  
  SUBROUTINE make_ham_x(mat1, id, N, outx)

    IMPLICIT NONE

    !Input: the Sigma-X Pauli's matrix and the identity
    TYPE(dmatrix), INTENT(IN) :: mat1, id

    !Output: the X-part of the Ising Model
    TYPE(dmatrix), INTENT(OUT) :: outx

    !Input: Number of particle
    INTEGER*8, INTENT(IN) :: N

    !Integer used in the DO cycles
    INTEGER*8 :: ii, jj, kk

    !Integer used to check the deallocation/allocation
    INTEGER*4 :: info

    !Temporary matrix used for the computation
    TYPE(dmatrix) :: trix1, trix2

    !I distinguish the case when N=2, i.e. two particles, and the case
    !in which N>2. The first case is simply and it consists in a tensor
    !product between the Sigma-X Pauli's matrices. The second case is
    !is similar to the Z-part and keep in consideration the three cases
    !in which I consider the first particle, then I consider the second
    !particle and finally the other.

    !Case in which N==2
    IF (N==2) THEN

       !Tensor product between Sigma-X and itself
       CALL t_prod(mat1, mat1, trix1)

       !Defining the dimension of the output matrix for N=2
       outx%row=2
       outx%col=2

       !Allocation of the output matrix
       ALLOCATE(outx%mat(2,2), STAT=info)

       !CHECKPOINT: checking the allocation
       CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it')

       !Defining the output matrix
       outx%mat=trix1%mat

       !Deallocation of the temporary matrix
       DEALLOCATE(trix1%mat, STAT=info)

       !CHECKPOINT: checking the deallocation
       CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it.')

       !Case in which N>2   
    ELSE

       !Inizialization of the output matrix's dimensions
       outx%row= INT((mat1%row)**N)
       outx%col=INT((mat1%col)**N)

       !Allocation of the ouput matrix
       ALLOCATE(outx%mat(outx%row, outx%col), STAT=info)

       !CHECKPOINT: checking the allocation
       CALL checkpoint(info/=0, 'Somenthing wrong in the allocation of the matrix&
            & in the X-part. Maybe not enough memory. Please, check it')

       !Initialization of the output matrix as null matrix
       outx%mat = (0d0, 0d0)

       !Computation of the X-part
       DO ii=1, N-1
          !CASE 1: first particle
          IF (ii==1) THEN

             !Tensor product between two Sigma-Xs matrix
             CALL t_prod(mat1, mat1, trix1)

             !Tensor products between the previous results and the
             !identity matrices.
             DO jj=1, N-2

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)

             ENDDO
             !CASE 2: second particle   
          ELSEIF (ii==2) THEN

             !Tensor product between the identity matrix and the Sigma-X
             CALL t_prod(id, mat1, trix1)

             !Tensor product between the previous result and the Sigma-X
             CALL t_prod(trix1, mat1, trix2)

             CALL substi(trix1, trix2)

             !Other tensor products
             DO jj=1, N-3

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)
             ENDDO
             !CASE 3: Other particles   
          ELSE

             !Initialization of the temporary matrix as 2X2 identity
             trix1%row=id%row
             trix1%col=id%col
             ALLOCATE(trix1%mat(trix1%row, trix1%col), STAT=info)

             !CHECKPOINT: checking the allocation
             CALL checkpoint(info/=0, 'Something wrong with the allocation of the&
                  & matrix. Please, check it')

             trix1%mat = id%mat

             !First ii-2 tensor products between the identities
             DO jj=1, ii-2

                CALL t_prod(trix1, id, trix2)

                CALL substi(trix1, trix2)

             ENDDO

             !Tensor product between the previous result and the Sigma-X
             CALL t_prod(trix1, mat1, trix2)

             CALL substi(trix1, trix2)

             !Tensor product between the previous result and the Sigma-X
             CALL t_prod(trix1, mat1, trix2)

             CALL substi(trix1, trix2)

             !If I am not considering the last particle, I must compute
             !the remeaning tensor products.
             IF (ii+1 /= N) THEN

                DO jj=1, N-ii-1

                   CALL t_prod(trix1, id, trix2)

                   CALL substi(trix1, trix2)

                ENDDO

             ENDIF


          ENDIF

          !Storing the output matrix
          DO jj=1, outx%row
             DO kk=1, outx%col
                outx%mat(jj,kk) = outx%mat(jj, kk) + trix1%mat(jj, kk)
             ENDDO
          ENDDO

          !Deallocation of the temporary matrix in order to restart
          !the cycle
          DEALLOCATE(trix1%mat, STAT=info)

          !CHECKPOINT: checking the deallocation
          CALL checkpoint(info/=0, 'Something wrogn with the deallocation. Please check it')

       ENDDO

    ENDIF

    !Computation of the trace of the output matrix
    CALL trace(outx)

    !Returning the result
    RETURN

  END SUBROUTINE make_ham_x
  !=====================================================================!
END MODULE ising
!======================================================================!


!======================================================================!
!This module contains the subroutines used for the RSRG algorithm to
!compute the approximated solution to the many-body problem.
!======================================================================!
MODULE reno

  !Using the previous modules
  USE tools
  USE ising

CONTAINS

  !====================================================================!
  !This subroutine computes the Identity matrix for the system with N
  !particles and N-1 particles. It takes as input the 2X2 identity
  !matrix and the size of the system. As output, it returns the
  !NxN identity and the N-1xN-1 identity.
  !====================================================================!
  SUBROUTINE identity(ide, N, identN, identNm1)

    !Input: 2x2 identity matrix
    TYPE(dmatrix), INTENT(IN) :: ide

    !Input: size of the system
    INTEGER*8, INTENT(IN) :: N

    !Output: NxN and N-1xN-1 identity matrices
    TYPE(dmatrix), INTENT(OUT) :: identN, identNm1

    !Integer used for the checkpoint and for the Do cycles
    INTEGER*8 :: info

    !Temporary matrix
    TYPE(dmatrix) :: trix1, trix2

    !CASE N=2
    IF(N==2)THEN

       !The N-1xN-1 identity matrix is simply the 2X2 identity matrix
       CALL copy(ide, identNm1)

       !The NxN identity matrix is the tensor product between the
       !2x2 identity matrices
       CALL t_prod(ide, identNm1, identN)
       !CASE N>2
    ELSE

       !Copying the 2x2 identity matrix into a temporary matrix
       CALL copy(ide, trix1)

       !Creation of the NxN identity 
       DO info=1, N-3

          !Tensor products
          CALL t_prod(trix1, ide, trix2)

          !Substitution of the temporary matrices
          CALL substi(trix1, trix2)

       ENDDO

       !Last tensor product for the N-1xN-1 identity matrix
       CALL t_prod(trix1, ide, identNm1)


       !Last two tensor products for the NxN identity matrix
       CALL t_prod(trix1, ide, trix2)

       !Substitution of the temporary matrices
       CALL substi(trix1, trix2)

       !Last tensor product
       CALL t_prod(trix1, ide, identN)

       !Deallocation of the temporary matrix
       DEALLOCATE(trix1%mat, STAT=info)

       !CHECKPOINT: checking the deallocation
       CALL checkpoint(info/=0, 'Somtheing wrong with the deallocation in ''idenitity''. Please, check it.')

    ENDIF

    !Returning the results
    RETURN

  END SUBROUTINE identity
  !====================================================================!


  !====================================================================!
  !This subroutine projects the matrix into the sub-system. As input,
  !it takes the matrix that will be projected, the projector. As output
  !it returns the projection of the matrix.
  !====================================================================!
  SUBROUTINE projection(matrice, proj)

    IMPLICIT NONE

    !Input and output: matrix that will be projected
    TYPE(dmatrix), INTENT(INOUT) :: matrice

    !Input: projector
    TYPE(dmatrix), INTENT(IN) :: proj

    !Temporary matrix
    TYPE(dmatrix) :: trix1, trix2

    !Integer used for the checkpoint
    INTEGER*8 :: info

    !Copying the input matrix into a temporary matrix
    CALL copy(matrice, trix1)

    !Deallocation of the input matrix
    DEALLOCATE(matrice%mat, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''projection''. &
         &Please, check it')

    !First projection
    CALL multi(trix1, proj, trix2)

    !Substitution of the temporary matrix
    CALL substi(trix1, trix2)


    !Copying the projector into the temporary matrix
    CALL copy(proj, trix2)

    !Compute the hermitian of the projector
    CALL hermi(trix2)    

    !Last multiplication
    CALL multi(trix2, trix1, matrice)

    !Deallocation of the temporary matrix
    DEALLOCATE(trix1%mat, trix2%mat, STAT=info)

    !CHECKPOINT: checking the deallocation
    CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''projection''. &
         &Please, check it.')

    !Returning the result
    RETURN

  END SUBROUTINE projection
  !====================================================================!

END MODULE reno
!======================================================================!



!======================================================================!  
!======================================================================!  
!======================================================================!  
!==========================MAIN PROGRAM================================!  
!======================================================================!  
!======================================================================!  
!======================================================================!  
PROGRAM renormalization

  USE tools
  USE ising
  USE reno

  IMPLICIT NONE

  !PAULI & IDENTITY
  TYPE(dmatrix) :: sigmax, sigmaz, id, identN, identNm1

  !Hamiltonian in the calculation
  TYPE(dmatrix) :: hamx, hamz, hamtot, hamint, &
       & ham2tot, h1, h2, h3, h4, h12, h23, h34, h2tot2, &
       & hamtotL, hamtotR, h1_2, h2_2, h12_2, h3_2, h4_2, h34_2

  !Interaction matrices
  TYPE(dmatrix) :: left, right

  !Temporary matrices
  TYPE(dmatrix) :: trix1, trix2, trix3, hamc

  !Density and reduced density matrices
  TYPE(dmatrix) :: rho1, rhoL, rhoR

  !Projector
  TYPE(dmatrix) :: projector

  !Discretization of the range of the strenght parameter and precision
  !used to evaluate the results
  REAL*8 :: dl, precision

  !Integer for the dimensions
  INTEGER*8 :: d, N, Nc, h

  !Integer for bins for the strenght parameter and for the range
  INTEGER*8 :: binlambda, range

  !Integer for the checkpoint
  INTEGER*8 :: info, limit

  !Integer for cycles
  INTEGER*8 :: ii, jj, kk, gg

  !Integer for the percentage and for the choose
  INTEGER*8 :: intero

  !Array used to store the ground states
  REAL*8, ALLOCATABLE :: levels(:)

  !Matrix to store the number of cycles and the number of particles
  REAL*8, ALLOCATABLE :: number(:,:), numbers(:,:)
  
  !Double precision for time spent
  REAL*8 :: ini, fine

  !Character for the choose
  CHARACTER :: cho, mano

  !Default parameters

  !Hilbert space dimension
  d=2
  !Size of the system
  N=2
  !Bins for the strenght parameter
  binlambda=200
  !Precision for the evaluation of the results
  precision=0.1
  !range for the variation of the strenght parameter
  range=3

  WRITE(*,'("=============================================================================")')
  WRITE(*,'("=============================================================================")')
  WRITE(*,'("=============================================================================")')
  WRITE(*, '("=====================Renormalization of the group flow======================")')
  WRITE(*,'("=============================================================================")')
  WRITE(*,'("=============================================================================")')
  WRITE(*,'("=============================================================================")')
  WRITE(*,*) ''
  WRITE(*,'("Do you want to proceed with the renormalization of the Ising Model? &
       &Type ''y'' to confirm or ''k'' to kill the program.")')
  WRITE(*,*)''
  READ*, mano
  WRITE(*,*)''
  

  IF(mano=='k')THEN
     STOP
  ELSEIF(mano=='y')THEN
     
     !The main program starts allowing to the user to change the default
     !Parameters
     WRITE(*, '("=========================== Renormalization Ising Model =====================")')
     WRITE(*,*)''
     WRITE(*,'("The default parameters are:")')
     WRITE(*,'("---------> N = ", I0)') N
     WRITE(*,'("---------> Binlambda = ", I0)') binlambda
     WRITE(*,'("---------> Precision = ", es20.10)') precision
     WRITE(*,'("---------> Lambda in [", F4.1,", ",F4.1, "]")') DBLE(0), DBLE(range)
     WRITE(*, *)''
     WRITE(*, '("Do you want change them? Type ''y'' for change or type ''k'' to kill the program.")')
     WRITE(*, *)''
     READ*, cho

     !If the user choose to change the default parameter, an infinite cycle
     !starts untill the user stops it.
     IF (cho=='y') THEN

        !Infinite cycle
        DO
           !Choose for which parameter change
           WRITE(*,*) ''
           WRITE(*, '("====================================================================")')
           WRITE(*,*) ''
           WRITE(*, '("What parameter do you wnat to change? Type the corrisponding number &
                &or type 0 to confirm the chosen or type -1 to exit from the program.")')
           WRITE(*,*)''
           WRITE(*,'("1) N = ", I0)') N
           WRITE(*,'("2) Binlambda = ", I0)') binlambda
           WRITE(*,'("3) Precision = ", es20.10)') precision
           WRITE(*,'("4) Lambda in [", F4.1,", ",F4.1, "]")') DBLE(0), DBLE(range)
           WRITE(*,*)''
           READ*, intero

           !If the user types '-1' the program stops itself
           IF (intero==-1) THEN
              STOP
              !If the user types '1', the program allows to change the size
              !of the system   
           ELSEIF (intero == 1) THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*, '("Please, insert the number of particles.")')
              WRITE(*,*)''
              READ*, N

              !Checking the input
              IF (N<2) THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a N bigger than 1. Setting N to the default parameter...")')
                 N=2
              ENDIF
              !If the user types '2', the programs allows to change the bins
              !for the strenght parameter
           ELSEIF(intero == 2) THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the number of bin for the parameter lambda.")')
              WRITE(*,*)''
              READ*, binlambda

              !checking the input
              IF (binlambda < 100) THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a bin for the parameter lambda bigger&
                      & than 100. Setting the bins to the default parameter...")')
                 binlambda=200
              ENDIF
              !If the user types '3', the programs allows to change the precision   
           ELSEIF(INTERO==3)THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the precision.")')
              WRITE(*,*)''
              READ*, precision

              !checking the input
              IF(precision > 1 .OR. precision< 0)THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a precision between 0 and 1. &
                      &Setting the precision to the default parameters...")')
                 precision = 0.1
              ENDIF
              !If the user types '4', the programs allows to change the range
              !of the strenght parameter   
           ELSEIF(intero==4)THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the total range for the parameter lambda.")')
              WRITE(*,*)''
              READ*, range
              !If the user types '0', the program starts to run
           ELSEIF(intero==0)THEN
              EXIT
           ELSE
              WRITE(*,*) ''
              WRITE(*,'("You make no choice.")')
           ENDIF
        ENDDO
     ELSEIF(cho == 'k') THEN
        STOP
     ENDIF

  ENDIF


  !Parameters used for the computation

  !Counter
  gg=0
  !Copy of the size of the system
  Nc=N
  !Integer used for the percentage
  intero=0
  !discretization of the strenght parameter
  dl = range/DBLE(binlambda-1)


  !Inizialization of the Pauli's matrices and the identity
  CALL pauli(sigmax, sigmaz, id, d)

  !Writing data on .gnu file
  IF(mano=='y')THEN
     OPEN(unit=100, file='gs-cycle_0.gnu')
     WRITE(100, '("set grid")')
     WRITE(100, '("set title ''GS/N as function of the cycles, with N=", I0, " and precision = ", es20.10, "''")')&
          &N, precision
     WRITE(100, '("set ylabel ''GS/N''")')
     WRITE(100, '("set xlabel ''# of cycles''")')
     WRITE(100, '("plot ''-'' u 1:2 w l lt 6 title ''GS/N with lambda = ", F3.1)')DBLE(0)
     
     OPEN(unit=200, file='gs-cycle_max.gnu')
     WRITE(200, '("set grid")')
     WRITE(200, '("set title ''GS/N as function of the cycles, with N=", I0, " and precision = ", es20.10, "''")')&
          &N, precision
     WRITE(200, '("set ylabel ''GS/N''")')
     WRITE(200, '("set xlabel ''# of cycles''")')
     WRITE(200, '("plot ''-'' u 1:2 w l lt 2 title ''GS/N with lambda = ", F4.1)')DBLE(range)
     
     OPEN(unit=300, file='l-gs.gnu')
     WRITE(300, '("set grid")')
     WRITE(300, '("set title ''GS/N of convergence as function of lambda&
          &, with N=2 and precision", es20.10, "''")')precision
     WRITE(300, '("set ylabel ''GS/N''")')
     WRITE(300, '("set xlabel ''lambda''")')
     WRITE(300, '("plot ''-'' u 1:2 w l lt 4 title ''GS/N''")')
     
     OPEN(unit=60, file='gs-c.gnu')
     WRITE(60, '("set grid")')
     WRITE(60, '("set title ''GS of convergence as function of # cycle &
          &, with N=2 and precision", es20.10, "''")')precision
     WRITE(60, '("set ylabel ''GS/N''")')
     WRITE(60, '("set xlabel ''cycles''")')
     WRITE(60, '("plot ''-'' u 1:2 w l lt 4 title ''GS''")')

  ENDIF
  
  !Time of the starting calculation
  CALL cpu_time(ini)  

  !Computation of the NxN identity and the N-1xN-1 identity
  CALL identity(id, N, identN, identNm1)

  
  !Computing the Z-part of the Ising model
  CALL make_ham_z(sigmaz, id, N, hamz)

  !Computing the X-part of the Ising model
  CALL make_ham_x(sigmax, id, N, hamx)

  !Computation as the strenght parameter varying
  DO ii=1, binlambda

     IF(mano=='y')THEN
        !IF statement for printing the  percentage
        IF ( INT((DBLE(ii)/DBLE(binlambda))*100d0) == intero) THEN
           WRITE(*, '("---> ",I0, "%")')intero
           intero = intero +1
        ENDIF
     ENDIF
     

     !Defining the dimension of the total Hamiltonian of the Ising Model
     hamtot%row=INT(d**N)
     hamtot%col=INT(d**N)

     !Allocation of the total Hamiltonian of the Ising Model
     ALLOCATE(hamtot%mat(hamtot%row, hamtot%col), STAT=info)

     !CHECKPOINT: checking the allocation
     CALL checkpoint(info/=0, 'Something wrong with the allocation the total &
          &Hamiltonian of the Ising Model. Please, check it.')

     !Defining the total Hamiltonian of the Ising Model
     hamtot%mat = (dl*(DBLE(ii-1)))*hamz%mat + hamx%mat

     !Computing the trace of the Total Hamiltonian
     CALL trace(hamtot)

     !==========================================================!
     !======================DOUBLING PART=======================!
     !==========================================================!

     !L-interaction term
     CALL t_prod(identNm1, sigmax, left)

     !computing the trace of the L-interaction term
     CALL trace(left)

     !R-interaction term
     CALL t_prod(sigmax, identNm1, right)

     !computing the trace of the R-interaction term
     CALL trace(right)

     !Allocation of the energies array
     ALLOCATE(levels(10000), STAT=info)

     !CHECKPOINT: checking the allocation
     CALL checkpoint(info/=0, 'Something wrong with the allocation of &
          &the array for the energies. Please, check it.')

     !Defining the default parameter

     !Counter to zero
     limit=0

     !Steps
     gg=0


     !Infinite loop for the evaluation of the ground state
     DO WHILE (limit==0)

        !Increasing the steps
        gg=gg+1

        !Defining the copy of the size and the double of itself
        Nc = N* (2**gg)

        !Tensor product between the L-interaction and the R-interaction
        !terms in order to gets the interaction Hamiltonian
        CALL t_prod(left, right, hamint)


        !Computin the trace of the Interaction term
        CALL trace(hamint)

        !Left-system
        CALL t_prod(hamtot, identN, trix1)

        !Right-system
        CALL t_prod(identN, hamtot, trix2)

        !Checking the correct dimension
        CALL checkpoint(SIZE(trix1%mat, 1) /= SIZE(trix2%mat, 1) .OR. &
             & SIZE(trix1%mat, 2) /= SIZE(trix2%mat,2) .OR. &
             & SIZE(trix1%mat, 1) /= SIZE(hamint%mat,1) .OR. &
             & SIZE(trix1%mat, 2) /= SIZE(hamint%mat,2), &
             &'The dimension of the three part are not equal. Please, check it.')

        !Defining the dimensions of the doubled Hamiltonian
        ham2tot%row = SIZE(hamint%mat,1)
        ham2tot%col = SIZE(hamint%mat,2)

        !Allocation of the doubled Hamiltonian
        ALLOCATE(ham2tot%mat(ham2tot%row, ham2tot%col), STAT=info)


        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong in the allocation of the &
             &doubled Hamiltonian. Please, check it.')

        !Defining the doubled hamiltonian
        ham2tot%mat = trix1%mat + trix2%mat + hamint%mat

        !Checking if the doubled Hamiltonian is hermitian
        DO jj=1, ham2tot%row
           DO kk=1, ham2tot%col
              CALL checkpoint( ABS(ham2tot%mat(jj,kk) - CONJG(ham2tot%mat(kk,jj))) > 0.1, &
                   &'The double hamiltonian is not hermitian. Please, check it')
           ENDDO
        ENDDO

        !Deallocation of the temporary matrix
        DEALLOCATE(trix1%mat, trix2%mat, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation in the main program.&
             & Please, check it.')

        !Copying the doubled Hamiltonian into a temporary matrix
        CALL copy(ham2tot, hamc)

        !Eigenproblem of the doubled Hamiltonian
        CALL eigenp(hamc)

        !Defining the E0 as the GS divided for the Number of particles
        levels(gg) = hamc%eigval(1)/Nc

        IF(mano=='y')THEN
           IF(ii==1)THEN
              WRITE(100, *) gg, levels(gg)
           ELSEIF(ii==binlambda)THEN
              WRITE(200,*) gg, levels(gg)
           ENDIF
        ENDIF

        !Checking the convergence
        IF(gg>1)THEN

           !If the previous E0 is similar to the actual, the program
           !stop the DO WHILE cycle
           IF(ABS(levels(gg-1)-levels(gg)) < precision) THEN

              WRITE(*,'("Convergence at step ", I0," for the cyle ", I0)')gg, ii
              !Stetting the counter to 1 in order to stop the DO WHILE
              !cycle
              limit=1

              IF(mano=='y')THEN
                 
                 !Writing the data in a .gnu file
                 WRITE(300,*) dl*(DBLE(ii-1)), levels(gg)
                 WRITE(60,*) gg, levels(gg)
                 
              ENDIF
              

              !If the convergence do not exist, the program stops the DO
              !WHILE cycle after 60 iteration
           ELSEIF (gg > 60) THEN

              !Counter to 1
              limit=1
           ENDIF
        ENDIF


        !Defining the dimension of the projectors
        projector%row=INT(2**(2*N))
        projector%col=INT(d**N)

        !ALLocation of the projector
        ALLOCATE(projector%mat(projector%row, projector%col), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation of the &
             &projector. Please, check it.')

        !Defining the projector
        DO jj=1, SIZE(projector%mat,1)
           DO kk=1, SIZE(projector%mat,2)
              projector%mat(jj, kk) = hamc%mat(jj, kk)
           ENDDO
        ENDDO

        !Projection of the doubled Hamiltonian
        CALL projection(ham2tot, projector)

        !substitution of the matrices: I substitute the projected
        !doubled Hamiltonian with the total Hamiltonian in order
        !to restart the computation
        CALL substi(hamtot, ham2tot)

        !Tensor product between the NxN identity and the L-interaction
        !part
        CALL t_prod(identN, left, trix1)

        !Substitution of the temporary matrix with the L-interaction
        !matix
        CALL substi(left, trix1)

        !Projection of the L-interaction part (of initial size 2N)
        CALL projection(left, projector)

        !Tensor product between the NxN identity and the R-interaction
        !part
        CALL t_prod(right, identN, trix1)

        !Substitution of the temporary matrix with the R-interaction
        !matix
        CALL substi(right, trix1)
        
        
        !Projection of the R-interaction part (of initial size 2N)
        CALL projection(right, projector)


        !Deallocation of the arrays used in the Do WHILE cycle
        DEALLOCATE(hamc%mat, hamc%eigval, projector%mat, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation of the matrices&
             & used in the DO WHILE cycle. Please, check it.')

     ENDDO

     !Default values

     !Size of the system to the starting one
     Nc=N

     !Deallocation of the matrices used in the DO cycles
     DEALLOCATE(hamtot%mat, left%mat, right%mat, levels, STAT=info)

     !CHECKPOINT: checking the deallocation
     CALL checkpoint(info/=0, 'Something wrong with the deallocation of the matrices&
          & used in the DO cycle. Please, check it.')

  ENDDO

  IF (mano=='y') THEN
     !Time for the end of computation
     CALL cpu_time(fine)

     !Defining the total time
     ini=fine-ini

     !If statement for the printing of the right time unit
     IF (ini<60.0) THEN
        PRINT*, 'Time:  ', ini, 'seconds'
     ELSEIF(ini>=60.0 .AND. ini<3600.0)THEN
        PRINT*, 'Time:  ', ini/60.0, 'minuts'
     ELSE
        PRINT*, 'Time:  ', ini/3600.0, 'hours'
     ENDIF
     
     CLOSE(100)
     CLOSE(200)
     CLOSE(300)
     CLOSE(60)

     CALL execute_command_line("gnuplot gs-cycle_0.gnu -p")

     CALL execute_command_line("gnuplot gs-cycle_max.gnu -p")
     
     CALL execute_command_line("gnuplot l-gs.gnu -p")
     
     CALL execute_command_line("gnuplot gs-c.gnu -p")
     
  ENDIF


 

  WRITE(*,'("Do you want to proceed with the Infinite DMRG? &
       &Type ''y'' to confirm or ''k'' to kill the program.")')
  WRITE(*,*)''
  READ*, mano
  WRITE(*,*)''


  IF(mano=='k')THEN
     STOP
  ELSEIF(mano=='y')THEN
     !The main program starts allowing to the user to change the default
     !Parameters
     WRITE(*, '("=========================== INFINITE DMRG =====================")')
     WRITE(*,*)''
     WRITE(*,'("The default parameters are:")')
     WRITE(*,'("---------> N = ", I0)') N
     WRITE(*,'("---------> Binlambda = ", I0)') binlambda
     WRITE(*,'("---------> Precision = ", es20.10)') precision
     WRITE(*,'("---------> Lambda in [", F4.1,", ",F4.1, "]")') DBLE(0), DBLE(range)
     WRITE(*, *)''
     WRITE(*, '("Do you want change them? Type ''y'' for change or type ''k'' to kill the program.")')
     WRITE(*, *)''
     READ*, cho

     !If the user choose to change the default parameter, an infinite cycle
     !starts untill the user stops it.
     IF (cho=='y') THEN

        !Infinite cycle
        DO
           !Choose for which parameter change
           WRITE(*,*) ''
           WRITE(*, '("====================================================================")')
           WRITE(*,*) ''
           WRITE(*, '("What parameter do you wnat to change? Type the corrisponding number &
                &or type 0 to confirm the chosen or type -1 to exit from the program.")')
           WRITE(*,*)''
           WRITE(*,'("1) N = ", I0)') N
           WRITE(*,'("2) Binlambda = ", I0)') binlambda
           WRITE(*,'("3) Precision = ", es20.10)') precision
           WRITE(*,'("4) Lambda in [", F4.1,", ",F4.1, "]")') DBLE(0), DBLE(range)
           WRITE(*,*)''
           READ*, intero

           !If the user types '-1' the program stops itself
           IF (intero==-1) THEN
              STOP
              !If the user types '1', the program allows to change the size
              !of the system   
           ELSEIF (intero == 1) THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*, '("Please, insert the number of particles.")')
              WRITE(*,*)''
              READ*, N

              !Checking the input
              IF (N<2) THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a N bigger than 1. Setting N to the default parameter...")')
                 N=2
              ENDIF
              !If the user types '2', the programs allows to change the bins
              !for the strenght parameter
           ELSEIF(intero == 2) THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the number of bin for the parameter lambda.")')
              WRITE(*,*)''
              READ*, binlambda

              !checking the input
              IF (binlambda < 100) THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a bin for the parameter lambda bigger&
                      & than 100. Setting the bins to the default parameter...")')
                 binlambda=200
              ENDIF
              !If the user types '3', the programs allows to change the precision   
           ELSEIF(INTERO==3)THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the precision.")')
              WRITE(*,*)''
              READ*, precision

              !checking the input
              IF(precision > 1 .OR. precision< 0)THEN
                 WRITE(*,*)''
                 WRITE(*,'("Wrong value. Insert a precision between 0 and 1. &
                      &Setting the precision to the default parameters...")')
                 precision = 0.01
              ENDIF
              !If the user types '4', the programs allows to change the range
              !of the strenght parameter   
           ELSEIF(intero==4)THEN
              WRITE(*,*) ''
              WRITE(*, '("====================================================================")')
              WRITE(*,*) ''
              WRITE(*,'("Please, insert the total range for the parameter lambda.")')
              WRITE(*,*)''
              READ*, range
              !If the user types '0', the program starts to run
           ELSEIF(intero==0)THEN
              EXIT
           ELSE
              WRITE(*,*) ''
              WRITE(*,'("You make no choice.")')
           ENDIF
        ENDDO
     ELSEIF(cho == 'k') THEN
        STOP
     ENDIF
  ENDIF
  
  !====================================================================!
  !====================================================================!
  !====================================================================!
  !===========================INFINITE DMRG============================!
  !====================================================================!
  !====================================================================!
  !====================================================================!
  IF(mano == 'y') THEN

     OPEN(unit=10, file='dmrg_gs-cycle_0.gnu')
     WRITE(10, '("set grid")')
     WRITE(10, '("set title ''GS/N as function of the cycles, with N=", I0, " and precision = ", es20.10, "''")')&
          &N, precision
     WRITE(10, '("set ylabel ''GS/N''")')
     WRITE(10, '("set xlabel ''# of cycles''")')
     WRITE(10, '("plot ''-'' u 1:2 w l lt 6 title ''GS/N with lambda = ", F3.1)')DBLE(0)

     OPEN(unit=20, file='dmrg_gs-cycle_max.gnu')
     WRITE(20, '("set grid")')
     WRITE(20, '("set title ''GS/N as function of the cycles, with N=", I0, " and precision = ", es20.10, "''")')&
          &N, precision
     WRITE(20, '("set ylabel ''GS/N''")')
     WRITE(20, '("set xlabel ''# of cycles''")')
     WRITE(20, '("plot ''-'' u 1:2 w l lt 2 title ''GS/N with lambda = ", F4.1)')DBLE(range)

     OPEN(unit=30, file='dmrg_l-gs.gnu')
     WRITE(30, '("set grid")')
     WRITE(30, '("set title ''GS/N of convergence as function of lambda&
          &, with N=2 and precision", es20.10, "''")')precision
     WRITE(30, '("set ylabel ''GS/N''")')
     WRITE(30, '("set xlabel ''lambda''")')
     WRITE(30, '("plot ''-'' u 1:2 w l lt 4 title ''GS/N''")')

     OPEN(unit=40, file='dmrg_c-l.gnu')
     WRITE(40, '("set grid")')
     WRITE(40, '("set title ''Number of cycles as function of lambda&
          &, with N=2 and precision", es20.10, "''")')precision
     WRITE(40, '("set ylabel ''# cycle''")')
     WRITE(40, '("set xlabel ''lambda''")')
     WRITE(40, '("plot ''-'' u 1:2 w l lt 1 title ''# of cycle''")')

     OPEN(unit=50, file='dmrg_gs-c.gnu')
     WRITE(50, '("set grid")')
     WRITE(50, '("set title ''GS of convergence as function of # cycle &
          &, with N=2 and precision", es20.10, "''")')precision
     WRITE(50, '("set ylabel ''GS/N''")')
     WRITE(50, '("set xlabel ''cycles''")')
     WRITE(50, '("plot ''-'' u 1:2 w l lt 4 title ''GS''")')

     !Allocation of the Total Hamiltonian of the ISing Model
     ALLOCATE(hamtot%mat(SIZE(hamx%mat, 1), SIZE(hamx%mat, 2)), STAT=info)

     !CHECKPOINT: checking the deallocation
     CALL checkpoint(info/=0, 'Something wrong with the allocation of Ham_tot in &
          &the infinte DMRG algorithm. Please, check it.')

     !Default values of the size
     Nc = N

     !Integer for the percentage
     intero=0

     CALL cpu_time(ini)
    
     !starting the calculation
     DO ii=1, binlambda
        
        !IF statement for printing the  percentage
        IF ( INT((DBLE(ii)/DBLE(binlambda))*100d0) == intero) THEN
           WRITE(*, '("---> ",I0, "%")')intero
           intero = intero +1
        ENDIF

        !Defining the left Hamiltonian as the Hamiltonian for the Ising Model
        hamtot%mat = dl*(DBLE(ii-1))*hamz%mat + hamx%mat !4x4

        !Copying the total Hamiltonian into temporary matrices
        CALL copy(hamtot, hamtotL)
        CALL copy(hamtot, hamtotR)     

        !Creation of the left-part of the interaction
        !LEFT_(2^N x 2^N) = 1_(2^(N-1) x 2^(N-1)) x S_x_(2x2)
        CALL t_prod(identNm1, sigmax, left)  !4x4

        !Creation of the right-part of the interaction
        !RIGHT_(2^N x 2^N) ) = S_x_(2x2) x 1_(2^(N-1) x 2^(N-1))
        CALL t_prod(sigmax, identNm1, right) !4x4

        !Allocation of the array to store the GS
        ALLOCATE(levels(1000000), STAT=info)

        !CHECKPOINT: checking the allocation.
        CALL checkpoint(info/=0, 'Something wrong with the allocation of the energies in &
             & the INFINITE DMRG part. Please, check it')

        !Setting the parameters before the DO WHILE cycle

        !Counter
        limit=0

        !Number of step
        gg=0

        !Initial size
        Nc=N

        !DO WHILE cycle untill the ''limit'' is equal to 0
        DO WHILE(limit /= 1)

           !Increasing the number of steps
           gg=gg+1

           !New size of the system
           Nc = 2*N + 2*gg

           !===============================================================!
           !======================NON-INTERACTION PART=====================!
           !===============================================================!

           !---------------------------------------------------------------!
           !Creating H1_(2^(2N+2) x 2^(2N+2)) =
           !   =  H_tot_(2^N x 2^N)  x  1_(2^(N+2) x 2^(N+2))

           !1) TEMP_(2^(N+1)  x  2^(N+1)) = H_tot_(2^N x 2^N)  x  1_(2x2)
           CALL t_prod(hamtotL, id, h1)

           !Copying the H1 into a temporary matrix in order to project it
           CALL copy(h1, h1_2) 

           !2) 1_(2^(N+1)  x  2^(N+1)) = 1_(2x2) x 1_(2^N x 2^N)
           CALL t_prod(id, identN, trix1) !8x8 
           
           !3) H1_(2^(2N+2)   x  2^(2N+2)) =
           !   =  TEMP_(2^(N+1)  x  2^(N+1))  x  1_(2^(N+1)  x  2^(N+1))
           CALL t_prod(h1, trix1, trix2)  !64x64

           !substitution of the temporary matrix with H1
           CALL substi(h1, trix2)
           

           !--------------------------------------------------------------!
           !Creating H2_(2^(2N+2)  x  2^(2N+2)) =
           !   = 1_(2^N x 2^N)  x  S_z_(2x2)  x 1_(2^(N+1) x 2^(N+1))

           !N.B. Beacuse there is only 1 particle, the interaction term in
           !the Ising Model do not appear. Thus, there is only
           !lambda*Sigma_z

           !1) TEMP_(2^(N+1) x 2^(N+1)) = 1_(2^N x 2^N)  x  S_z_(2x2)
           CALL t_prod(identN, sigmaz, trix2) !8x8

           !1.1) lambda*TEMP_(2^(N+1) x 2^(N+1))
           trix2%mat = dl*(DBLE(ii-1))*trix2%mat

           !Copying the temporary matrix into another temporary matrix in
           !order to project it
           CALL copy(trix2, h2_2)

           !2) H2_(2^(2N+2) x 2^(2N+2)) =
           !    =  lambda*TEMP_(2^(N+1) x 2^(N+1)) x 1_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix2, trix1, h2) !64x64

           !CHECKPOINT: Checking the correct dimension of the Hamiltonans
           CALL checkpoint(SIZE(h1%mat, 1) /= SIZE(h2%mat, 1), &
                &'The number of rows of H1 and H2 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h1%mat, 2) /= SIZE(h2%mat, 2), &
                &'The number of columns of H1 and H2 are not equal. Plase, check it.')

           !Defining the dimension of the H_(2^(2N+2) x 2^(2N+2)) for the
           !infinite DMRG
           h2tot2%row=SIZE(h1%mat,1)
           h2tot2%col=SIZE(h1%mat,2)

           !Allocation of the Hamiltonian
           ALLOCATE(h2tot2%mat(h2tot2%row, h2tot2%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong with the allocation of the H_(2N+2). Please, check it.')

           !Defining the Hamiltonian for the Infinite DMRG
           h2tot2%mat = h1%mat + h2%mat

           !DEALLOCATION of the temporary matrix
           DEALLOCATE(trix2%mat, h1%mat, h2%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H2''. Please, check it.')


           !--------------------------------------------------------------!
           !Creating H3_(2^(2N+2)  x  2^(2N+2)) =
           !   = 1_(2^(N+1) x 2^(N+1))  x  S_z_(2x2)  x 1_(2^N x 2^N)

           !N.B. Beacuse there is only 1 particle, the interaction term in
           !the Ising Model do not appear. Thus, there is only
           !lambda*Sigma_z

           !1) TEMP_(2^(N+1) x 2^(N+1)) = S_z_(2x2) x  1_(2^N x 2^N)
           CALL t_prod(sigmaz, identN, trix2) !8x8

           !1.1) lambda*TEMP_(2^(N+1) x 2^(N+1))
           trix2%mat = dl*(DBLE(ii-1))*trix2%mat

           !Copying the temporary matrix into another temporary matrix in
           !order to project it
           CALL copy(trix2, h3_2)

           !2) H3_(2^(2N+2) x 2^(2N+2)) =
           !    =  1_(2^(N+1) x 2^(N+1)) x lambda*TEMP_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix1, trix2, h3) !64x64

           !Checking the dimension of the previo
           CALL checkpoint(SIZE(h2tot2%mat, 1) /= SIZE(h3%mat, 1), &
                &'The number of rows of H2 and H3 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h2tot2%mat, 2) /= SIZE(h3%mat, 2), &
                &'The number of columns of H2 and H3 are not equal. Plase, check it.')

           !Summing the new Hamiltonian
           h2tot2%mat = h2tot2%mat + h3%mat

           !DEALLOCATION of the temporary matrix
           DEALLOCATE(trix2%mat, h3%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H3''. Please, check it.')

           !--------------------------------------------------------------!
           !Creating H4_(2^(2N+2) x 2^(2N+2)) =
           !    = 1_(2^(N+1) x 2^(N+1) x 1_(2x2) x H_tot_(2^N x 2^N)

           !1) TEMP_(2^(N+1) x 2^(N+1)) = 1_(2x2) x H_tot_(2^N x 2^N)
           CALL t_prod(id, hamtotR, trix2) !8x8

           !Copying the temporary matrix into another temporary matrix in
           !order to project it
           CALL copy(trix2, h4_2)

           !2) H4_(2^(2N+2) x 2^(2N+2)) =
           !    = 1_(2^(N+1) x 2^(N+1)) x TEMP_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix1, trix2, h4) !64x64

           !Checking the dimension of the Hamiltonians
           CALL checkpoint(SIZE(h2tot2%mat, 1) /= SIZE(h4%mat, 1), &
                &'The number of rows of H3 and H4 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h2tot2%mat, 2) /= SIZE(h4%mat, 2), &
                &'The number of columns of H2 and H4 are not equal. Plase, check it.')

           !Summing the new hamiltonian
           h2tot2%mat = h2tot2%mat + h4%mat

           !DEALLOCATION of the temporary matrix
           DEALLOCATE(trix2%mat, h4%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H4''. Please, check it.')       


           !===============================================================!
           !========================INTERACTION PART=======================!
           !===============================================================!

           !--------------------------------------------------------------!
           !Creating H_int_12_(2^(2N+2) x 2^(2N+2)) =
           !     = 1_(2^(N-1) x 2^(N-1)) x S_x_(2x2) x S_x_(2x2) x
           !                     x  1_(2^(N+1) x 2^(N+1))

           !1) TEMP_(2^(N+1) x 2^(N+1)) = LEFT_(2^N x 2^N) x S_x_(2x2)
           CALL t_prod(left, sigmax, trix2) !8x8

           !copying the temporary matrix into another temporary matrix in
           !order to project it
           CALL copy(trix2, h12_2)

           !2) H_int_12_(2^(2N+2) x 2^(2N+2)) =
           !     = TEMP_(2^(N+1) x 2^(N+1)) x 1_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix2, trix1, h12)   !64x64

           !Checking the dimension of the previous Hamiltonian
           CALL checkpoint(SIZE(h2tot2%mat, 1) /= SIZE(h12%mat, 1), &
                &'The number of rows of H4 and H12 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h2tot2%mat, 2) /= SIZE(h12%mat, 2), &
                &'The number of columns of H4 and H12 are not equal. Plase, check it.')

           !Summing the new Hamiltonian
           h2tot2%mat = h2tot2%mat + h12%mat

           !Deallocation of the temporary matrix
           DEALLOCATE(trix2%mat, h12%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H_12''. Please, check it.')

           !--------------------------------------------------------------!
           !Creating H_int_23_(2^(2N+2) x 2^(2N+2)) =
           !     = 1_(2^(N) x 2^(N)) x S_x_(2x2) x S_x_(2x2) x
           !                     x  1_(2^(N) x 2^(N))

           !1) TEMP_1_(2^(N+1) x 2^(N+1)) = 1_(2^N x 2^N) x S_x_(2x2)
           CALL t_prod(identN, sigmax, trix2) !8x8

           !2) TEMP_2_(2^(N+1) x 2^(N+1)) = S_x_(2x2) x 1_(2^N x 2^N)
           CALL t_prod(sigmax, identN, trix3) !8x8

           !3) H_int_23_(2^(2N+2) x 2^(2N+2)) =
           !      = TEMP_1_(2^(N+1) x 2^(N+1)) x TEMP_2_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix2, trix3, h23) !64x64

           !Checking the dimension of the previous Hamiltonian
           CALL checkpoint(SIZE(h2tot2%mat, 1) /= SIZE(h23%mat, 1), &
                &'The number of rows of H12 and H23 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h2tot2%mat, 2) /= SIZE(h23%mat, 2), &
                &'The number of columns of H12 and H23 are not equal. Plase, check it.')

           !Summing the new Hamiltonian
           h2tot2%mat = h2tot2%mat + h23%mat

           !Deallocation of the temporary matrices
           DEALLOCATE(trix2%mat, trix3%mat, h23%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H_23''. Please, check it.')
           !--------------------------------------------------------------!
           !Creating H_int_34_(2^(2N+2) x 2^(2N+2)) =
           !     = 1_(2^(N+1) x 2^(N+1)) x S_x_(2x2) x S_x_(2x2) x
           !                     x  1_(2^(N-1) x 2^(N-1))

           !1) TEMP_(2^(N+1) x 2^(N+1)) = S_x_(2x2) x RIGHT_(2^N x 2^N)
           CALL t_prod(sigmax, right, trix2) !8x8

           !Copying the temporary matrix into another temporary matrix in
           !order to project it
           CALL copy(trix2, h34_2)

           !2) H_int_34_(2^(2N+2) x 2^(2N+2)) =
           !    = 1_(2^(N+1) x 2^(N+1)) x TEMP_(2^(N+1) x 2^(N+1))
           CALL t_prod(trix1, trix2, h34)   !64x64

           !Checking the dimension of the previous Hamiltonian
           CALL checkpoint(SIZE(h23%mat, 1) /= SIZE(h34%mat, 1), &
                &'The number of rows of H23 and H34 are not equal. Plase, check it.')
           CALL checkpoint(SIZE(h23%mat, 2) /= SIZE(h34%mat, 2), &
                &'The number of columns of H2 and H4 are not equal. Plase, check it.')

           h2tot2%mat = h2tot2%mat + h34%mat

           !Deallocation of the temporary matrix
           DEALLOCATE(trix2%mat, h34%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation in ''H_34''. Please, check it.')


           !==============================================================!
           !===============DIAGONALIZATION OF H_(2N+2)====================!
           !==============================================================!

           !Allocation of the array where the eigenvalues will be stored
           ALLOCATE(h2tot2%eigval(h2tot2%row), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong with the allocation of the array &
                for the eigenvalues in the Infinite DMRG. Please, chek it.')

           !Checking if the H_(2N+2) Hamiltonian is hermitian
           DO jj=1, h2tot2%row
              DO kk=1, h2tot2%col
                 CALL checkpoint( ABS(h2tot2%mat(jj,kk) - CONJG(h2tot2%mat(kk,jj))) > 0.1, &
                      &'The H_(2N+2) Hamiltonian is not hermitian. Please, check it')
              ENDDO
           ENDDO

           !Eigenproblem for the H_(2N+2) Hamiltonian
           CALL eigenp(h2tot2)

           !Defining the GS as the total GS dived by the Number of particles
           levels(gg) = h2tot2%eigval(1)/Nc
           
           IF(ii==1)THEN
              WRITE(10, *) gg, levels(gg)
           ELSEIF(ii==binlambda)THEN
              WRITE(20,*) gg, levels(gg)
           ENDIF
           
           !Checking the convergence
           IF(gg>1)THEN
              IF(ABS(levels(gg) - levels(gg-1)) < precision)THEN
                 WRITE(*,'("Convergence at step ", I0," for the cyle ", I0)')gg, ii
                 WRITE(30,*) dl*(DBLE(ii-1)), levels(gg)
                 WRITE(40,*) dl*(DBLE(ii-1)), gg
                 WRITE(50,*) gg, levels(gg)
                 limit=1
              ENDIF
           ENDIF

           !----------------------Density matrix for the GS----------------!

           !Temporary matrices to make the first eigenvector a matrix
           !The ket is a 2N+2 x 1 matrix
           trix2%row = h2tot2%row
           trix2%col = 1

           !The bra is a 1x 2N+2 matrix
           trix3%row = 1
           trix3%col = h2tot2%col

           !Allocation of the bra and ket
           ALLOCATE(trix2%mat(trix2%row, trix2%col), trix3%mat(trix3%row, trix3%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong with the allocation of the bra and &
                &the ket. Please, check it.')

           !Defining the ket
           DO jj=1, trix2%row
              trix2%mat(jj,1) = h2tot2%mat(jj,1)
           ENDDO

           !Defining the bra
           DO jj=1, trix3%col
              trix3%mat(1, jj) = h2tot2%mat(jj, 1)
           ENDDO

           !Tensor product between the bra and the ket in order to get the
           !density matrix
           CALL t_prod(trix3, trix2, rho1)

           !Deallocation of the matrices used in the DO WHILE cycle
           DEALLOCATE(h2tot2%mat, h2tot2%eigval, trix2%mat, trix1%mat, trix3%mat,&
                & hamtotL%mat, hamtotR%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong with the deallocation of the &
                &matrices used in the DO WHILE cycle in DMRG. Please, check it.')

           !Checking if the trace of Rho is 1
           CALL trace(rho1)

           !CHECKPOINT: checking the trace of the density matrix
           CALL checkpoint(ABS(ABS(rho1%trace))-1d0  >  0.001, 'The trace of &
                &the density matrix is not equal 1. Please, check it.')

           !CHECKPOINT: checking if the density matrix is hermitian
           DO jj=1, rho1%row
              DO kk=1, rho1%col
                 CALL checkpoint( ABS(rho1%mat(jj,kk) - CONJG(rho1%mat(kk,jj))) > 0.1, &
                      &'The density matrix is not hermitian. Please, check it')
              ENDDO
           ENDDO

           !---------------------REDUCED DENSITY MATRIX-------------------!

           !--------------------------------L-----------------------------!

           !Computing the reduced density matrix for the L-subsystem
           CALL reduced(rho1, rhoL, 1_8, 8_8)

           !Diagonalization of the reduced density matrix for L-sub-system
           CALL eigenp(rhoL)

           !Checking the eigenvalues:
           DO jj=1, SIZE(rhoL%eigval)

              !1) the small eigenvalues near to 0, will be set to 0
              IF(ABS(rhoL%eigval(jj) ) < 1E-8) THEN
                 rhoL%eigval(jj) = +0d0
              ENDIF

              !2) If the eigenvalues are positive
              CALL checkpoint(rhoL%eigval(jj) < 0d0, 'The eigenvalues of the reduced &
                   &density matrix for L are not positive. Please, check it.')
           ENDDO

           !3) the sum of the eigenvalues must be 1
           CALL checkpoint(ABS(SUM(rhoL%eigval) - 1d0) > 0.0001, &
                &'The sum of the eigenvalues for the L-sub-system is not 1. Please, check it.')

           !Building the projector
           projector%row=SIZE(rhoL%mat,1)
           projector%col=(d**N)

           !Allocation of the projector
           ALLOCATE(projector%mat(projector%row, projector%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong in the allocation of the projector&
                & for L in the DMRG. Please, check it.')

           !Building the projector in the descending order
           DO jj=1, projector%row
              DO kk=1, projector%col
                 projector%mat(jj,kk) = rhoL%mat(jj, (d**(N+1))-kk+1)
              ENDDO
           ENDDO

           !Defining the projected Hamiltonian for theIsing Model
           hamtotL%row=h1_2%row
           hamtotL%col=h1_2%col


           !Allocation of the Hamiltonian for the Ising Model
           ALLOCATE(hamtotL%mat(hamtotL%row, hamtotL%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong with the allocation of the &
                &hamtotL in DMRG. Please, check it.')

           !Defining the L-part of the total Hamiltonian
           hamtotL%mat = h1_2%mat + h2_2%mat + h12_2%mat

           !Projecting the L-part of the total Hamiltonian
           CALL projection(hamtotL, projector)

           !Tensor product between the 1_(2^N x 2^N) identity ad the S_x
           CALL t_prod(identN, sigmax, trix1)

           !Substituting the temporary matrix with the LEFT matrix
           CALL substi(left, trix1)

           !Projecting the LEFT matrix
           CALL projection(left, projector)


           !--------------------------------R-----------------------------!

           !Deallocation of the projector
           DEALLOCATE(projector%mat, STAT=info)

           !CHECKPOINT: checking the deallocation
           CALL checkpoint(info/=0, 'Something wrong in the deallocation of the &
                &projector in DMRG. Please, check it.')

           !Computing the reduced density matrix for the R-subsystem
           CALL reduced(rho1, rhoR, 2_8, 8_8)

           !Diagonalization of the reduced density matrix for R-sub-system
           CALL eigenp(rhoR)

           !Checking the eigenvalues:
           DO jj=1, SIZE(rhoR%eigval)

              !1) the small eigenvalues near to 0, will be set to 0
              IF(ABS(rhoR%eigval(jj) ) < 1E-8) THEN
                 rhoR%eigval(jj) = +0d0
              ENDIF

              !2) If the eigenvalues are positive
              CALL checkpoint(rhoR%eigval(jj) < 0d0, 'The eigenvalues of the reduced &
                   &density matrix for R are not positive. Please, check it.')
           ENDDO

           !3) the sum of the eigenvalues must be 1
           CALL checkpoint(ABS(SUM(rhoR%eigval) - 1d0) > 0.0001, &
                &'The sum of the eigenvalues for the R-sub-system is not 1. Please, check it.')

           !Building the projector
           projector%row=SIZE(rhoR%mat,1)
           projector%col=(d**N)

           !Allocation of the projector
           ALLOCATE(projector%mat(projector%row, projector%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong in the allocation of the projector&
                & for the R-part in the DMRG. Please, check it.')

           !Building the projector in the descending order
           DO jj=1, projector%row
              DO kk=1, projector%col
                 projector%mat(jj,kk) = rhoR%mat(jj, (d**(N+1))-kk+1)
              ENDDO
           ENDDO

           !Defining the projected Hamiltonian for theIsing Model
           hamtotR%row=h3_2%row
           hamtotR%col=h3_2%col

           !Allocation of the Hamiltonian for the Ising Model
           ALLOCATE(hamtotR%mat(hamtotR%row, hamtotR%col), STAT=info)

           !CHECKPOINT: checking the allocation
           CALL checkpoint(info/=0, 'Something wrong with the allocation of the &
                &hamtotR in DMRG. Please, check it.')

           !Defining the L-part of the total Hamiltonian
           hamtotR%mat = h3_2%mat + h4_2%mat + h34_2%mat

           !Projecting the L-part of the total Hamiltonian
           CALL projection(hamtotR, projector)

           !Tensor product between the S_x_(2x2) and the 1_(2^N x 2^N) identity
           CALL t_prod(sigmax, identN, trix1)

           !Substituting the temporary matrix with the RIGHT matrix
           CALL substi(right, trix1)

           !Projecting the LEFT matrix
           CALL projection(right, projector)

           DEALLOCATE(projector%mat)
           DEALLOCATE(rho1%mat, rhoL%mat, rhoL%eigval, rhoR%mat, rhoR%eigval,&
                &h1_2%mat, h2_2%mat, h3_2%mat, h4_2%mat, h12_2%mat, h34_2%mat, STAT=info)



        ENDDO

        !Deallocation of the matrices in the DO cycle
        DEALLOCATE(levels, hamtot%mat, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation of the &
             &matrices used in the DO cycle in DMRG. Please, check it.')

     ENDDO
     
     !Time for the end of computation
     CALL cpu_time(fine)

     !Defining the total time
     ini=fine-ini

     !If statement for the printing of the right time unit
     IF (ini<60.0) THEN
        PRINT*, 'Time:  ', ini, 'seconds'
     ELSEIF(ini>=60.0 .AND. ini<3600.0)THEN
        PRINT*, 'Time:  ', ini/60.0, 'minuts'
     ELSE
        PRINT*, 'Time:  ', ini/3600.0, 'hours'
     ENDIF


     CLOSE(10)
     CLOSE(20)
     CLOSE(30)
     CLOSE(40)
     CLOSE(50)
     
     CALL execute_command_line("gnuplot dmrg_gs-cycle_0.gnu -p")

     CALL execute_command_line("gnuplot dmrg_gs-cycle_max.gnu -p")
     
     CALL execute_command_line("gnuplot dmrg_l-gs.gnu -p")

     CALL execute_command_line("gnuplot dmrg_c-l.gnu -p")

     CALL execute_command_line("gnuplot dmrg_gs-c.gnu -p")

  ENDIF
  

  !Deallocation of the remaining matrices
  DEALLOCATE(sigmax%mat, sigmaz%mat, id%mat, identN%mat, identNm1%mat, STAT=info)

  !CHECKPOINT: checking the deallocation of the remaning matrices
  CALL checkpoint(info/=0, 'Something wrong with the deallocation of the main matrices. &
       &Please, check it.')  

  
  !Deallocation of the X- and Z-part of the Ising Model
  DEALLOCATE(hamz%mat, hamx%mat, STAT=info)

  !CHECKPOINT:  checking the deallocation
  CALL checkpoint(info/=0, 'Something wrong with the deallocation of the &
       &X- and Z-part of the Ising Model. Plase, check it.')
END PROGRAM renormalization
