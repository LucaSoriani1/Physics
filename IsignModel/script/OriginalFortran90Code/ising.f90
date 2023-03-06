!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!The program are made to compute the Hamiltonia operator for ''N''
!particles with spin 1/2 that form the Ising Model. In this model
!appears the Pauli's matrices and the interaction strength.
!The program computes the 2^N X 2^N representation for the Hamiltonian
!for different ''N'' and for the parameters between 0 and 3.
!Thus, the program compute the eigenvalues and eigenvectors problem.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!



!**********************************************************************!
!Main module that contains all the subroutine and the personalized
!variable type used for the computation.
MODULE ising

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
    CALL checkpoint(info/=0, 'Somenthing wrong with the deallocation. Please check it')

    !Defining the dimensions of the new matrix as the dimensions of
    !the second matrix
    mat1%row=mat2%row
    mat1%col=mat2%col

    !Allocation of the new matrix with the new dimension
    ALLOCATE(mat1%mat(mat1%row, mat1%col), STAT=info)

    !CHECKPOINT: checking the allocation of the copied matrix
    CALL checkpoint(info/=0, 'Somenthing wrong with the allocation. Please check it')

    !Copying the second matrix in the first
    mat1%mat = mat2%mat

    !Deallocation of the second matrix
    DEALLOCATE(mat2%mat, STAT=info)

    !CHECKPOINT: checking the deallocation of the copied matrix
    CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it')

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

    !CHECKPOINT: checking the allocation of the matrix
    CALL checkpoint(info/=0, 'Something wrong with the allocation. Please check it1')
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
    
  END SUBROUTINE t_prod
  !=====================================================================!
  

  
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

       outx%row=2
       outx%col=2

       ALLOCATE(outx%mat(2,2))

       outx%mat=trix1%mat

       DEALLOCATE(trix1%mat)

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
    
  END SUBROUTINE make_ham_x
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

  END SUBROUTINE trace
  !=====================================================================!
    
 
END MODULE ising


!=====================================================================!  
!=====================================================================!  
!==========================MAIN PROGRAM===============================!  
!=====================================================================!  
!=====================================================================!  
!=====================================================================!  
PROGRAM hami

  !Using the module
  USE ising

  IMPLICIT NONE

  !The types 'dmatrix' used in the computation
  TYPE(dmatrix) :: sigmax, sigmaz, ide, hamx, hamz, hamtot

  !Integer for the computation
  INTEGER*8 :: d, binlambda, N

  !Integer used in the DO cycles
  INTEGER*8 :: ii, jj, ll, kk

  !Double precision used for the discretization of lambda
  REAL*8 :: dl

  !Double precision used to computation of the time spent
  REAL*8 :: start, stopp, time1, time2

  !Eigenproblem variables
  INTEGER*4 :: lwork=-1, lwmax=1000000, info

  DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: work

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rwork

  !Character for the choice
  CHARACTER :: cho

  !Integer for the percentage
  INTEGER*4 :: inte

  !Integer for the first k level
  INTEGER*8 :: kappa
  
  !Matrix where store the eigenvalues
  REAL*8, ALLOCATABLE :: eigenvalues(:,:)

  !Setting the default parameters
  d=2
  binlambda=300
  N=3

  !N=10 circa 46 minuti
  !N=11 circa 4.30 ore
  !N=12 circa 12 ore
  !Allowing to the user to change the default parameters
  WRITE(*, '("===========================================================")')
  WRITE(*, *) '' 
  WRITE(*, '("The default parameter is:")')
  WRITE(*, '("N = "I0)') N
  WRITE(*, '("Do you want to change it? Type ''y'' to change or type ''k'' to stop the program")')
  WRITE(*,*)''
  READ*, cho

  !If the answer is 'y'. Changing the number of particle. The precondition is N>=1
  IF (cho == 'y') THEN
     WRITE(*,*) ''
     WRITE(*, '("===========================================================")')
     WRITE(*,*) ''
     WRITE(*, '(" -----> N = ", I0)') N
     WRITE(*,*) ''
     WRITE(*, '("Please, insert the number of particles")')
     WRITE(*,*) ''
     READ*, N


     !CHECKPOINT: checking the correct values of N
     CALL checkpoint(N<1, 'Please, insert a number of particle bigger then 0')

  ELSEIF(cho=='k')THEN
     STOP
  ENDIF

  !Percentage to 0
  inte=0

  !Discretization of the strength parameter
  dl = 3d0/DBLE(binlambda-1)

  !Defualt k-leves if N=1, N=2 ore N=3. If N>3, the user can chose
  !how many levels he wants to plot
  !If N=1, there are only a 2X2 matrix, thus 2 eigenvalues
  IF(N==1) THEN
     kappa=2
  !If N=2, there is a 4X4 matrix, thus 4 eigenvalues
  ELSEIF(N==2) THEN
     kappa=4
  !If N=8, there is a 8X8 matrix, thus 8 eigenvalues
  ELSEIF(N==3)THEN
     kappa=8
  !Else N>3 there is a 2^N X 2^N matrix, thus 2^N eigenvalues.
  !The infinite loop is used to check the correct value of k.   
  ELSE   
     !Infinite loop
     DO
        !Defining the first k-level of eigenvalues
        WRITE(*,*) ''
        WRITE(*, '("===========================================================")')
        WRITE(*,*) ''
        WRITE(*,'("How many levels do you want to print? (Between 2 and ", I0,")")')INT(d**N)
        WRITE(*,*)''
        READ*, kappa
        WRITE(*,*)''
        
        !Checking the kappa
        IF(kappa<2 .OR. kappa > INT(d**N))THEN
           WRITE(*, '("Wrong value. Please, insert another one between 2 and ", I0)') INT(d**N)
        ELSE
           EXIT
        ENDIF
     ENDDO
  ENDIF
  
  !Allocation of the matrix in which storing the eigenvalues
  ALLOCATE(eigenvalues(binlambda, kappa), STAT=info)

  !CHECKPOINT: checking the allocation
  CALL checkpoint (info/=0, 'Something wrong with the allocation. Please, check it.')

  !Opening the file in which store the data
  OPEN(unit=1, file='eigval.gnu')

  !Writing options in the file .gnu
  !If the user wants to plot more than 8 eigenvalues,
  !the program keeps clean the graph not labelling the plots
  IF (kappa>8) THEN
     WRITE(1, '("set nokey")')
  ENDIF
    WRITE(1, '("set grid")')
  WRITE(1, '("set title ''Evolution of the first four eigenvalues with N=", I0)')N
  WRITE(1, '("set ylabel ''Eigenvalue''")')
  WRITE(1, '("set xlabel ''Parameters lambda''")')
  WRITE(1, '("plot ''-'' u 1:2 w l title ''Eigenvalue k=1'', "A1)')CHAR(92)
  DO ii=1, kappa-2
     WRITE(1, '(" ''-'' u 1:2 w l title ''Eigenvalue k=",I0,"'', "A1)')ii+1, CHAR(92)
  ENDDO
  WRITE(1, '(" ''-'' u 1:2 w l title ''Eigenvalue k=",I0)')kappa


  !Starting time for the computation
  CALL cpu_time(start)

  !Initialization of the Pauli's matrix
  CALL pauli(sigmax, sigmaz, ide, d)


  !Starting of the cycles
  DO ll=1, binlambda

     !IF statement for the percentage
     IF ( INT((DBLE(ll)/DBLE(binlambda))*100d0) == inte) THEN
        WRITE(*, '("---> ",I0, "%")')inte
        inte = inte +1
      ENDIF

     !Case in which N=1, thus the interaction part is null.
     IF (N==1) THEN

        !If N=1, the Hamiltonian operator is a 2X2 matrix and
        !correspond to the Sigma-Z Pauli's matrix
        hamtot%row=d
        hamtot%col=d

        !Allocation of the total Hamiltonian
        ALLOCATE(hamtot%mat(d,d), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it.')

        !Defining the total hamiltonian as the Sigma-Z multiply the
        !lambda parameter
        hamtot%mat=dl*(ll-1)*sigmaz%mat

        !Allocation of the array for the eigenvalue and for the
        !computation of the eigenvalues.
        ALLOCATE(hamtot%eigval(d), work(lwmax), rwork(3*d-2), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it.')

        !Defining lwrok=-1 in order to compute the optimal 'lwork'
        !for the computation of the eigenvalues
        lwork=-1

        !Computing the optimal 'lwork' for the computation of the
        !eigenvalues
        CALL ZHEEV('N', 'U', hamtot%col, hamtot%mat, hamtot%row, hamtot%eigval, work, &
             & lwork, rwork, info)

        !CHECKPOINT: checking the working of the ZHEEV subroutine
        CALL checkpoint(info/=0, 'Something wrong with the first ZHEEV. Please, check it.')

        !Defining of the optimal 'lwork'
        lwork = MIN(lwmax, INT(work(1)))

        !Deallocation of the 'work' array
        DEALLOCATE(work, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please check it.')

        !Allocation of the 'work' array of the right size
        ALLOCATE(work(MAX(1, lwork)), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it.')

        !Computing of the eigenvalues
        CALL ZHEEV('N', 'U', hamtot%col, hamtot%mat, hamtot%row, hamtot%eigval, work, &
             & lwork, rwork, info)

        !CHECKPOINT: checking the work of the ZHEEV subroutine
        CALL checkpoint(info/=0, 'Something wrong with the second ZHEEV. Please, check it.')

        !Deallocation of the arrays used in the ZHEEV parameter
        DEALLOCATE(work, rwork, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it.')

        !Case in which N>1
     ELSE

        !In order to determine the maximum N for the computation, I
        !compute the time spent for the computation of the creation
        !of the Z-part. If this is bigger than 40s, the program
        !stops itself

        !Starting time for the computation of the Z-part
        CALL cpu_time(time1)

        !Computing the Z-part of the Hamiltonian operator
        CALL make_ham_z(sigmaz, ide, N, hamz)

        !Computing the trace of the Z-part matrix.
        CALL trace(hamz)


        !Ending time for the computation of the Z-part           
        CALL cpu_time(time2)

        !Calculation of the time spent for the computation
        time1=time2-time1

        !CHECKPOINT: if the time spent is bigger than 40 the program stops.
        !CALL checkpoint( time1>40.0, 'Computation time bigger than 40s. Too much particles.')

        !Computing the X-part of the Hamiltonian operator
        CALL make_ham_x(sigmax, ide, N, hamx)

        !Computing the trace of the X-part matrix.
        CALL trace(hamx)

        !Degining the dimension of the Hamiltonian operator
        hamtot%row=INT(d**N)
        hamtot%col=INT(d**N)

        !Allocation of the Hamiltonian operator
        ALLOCATE(hamtot%mat(hamtot%row, hamtot%col), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation of the Hamiltonian&
             & operator. Maybe not enough memory. Please, check it.')

        !Initialization of the Hamiltonian operator to null-matrix
        hamtot%mat = (0d0, 0d0)

        !Storing the Hamiltonian operator
        hamtot%mat = dl*(DBLE(ll-1))*hamz%mat + hamx%mat

        !Computation the trace of the Hamiltonian operator
        CALL trace(hamtot)

        !Allocation of the eigenvalues array and the arrays used in
        !the ZHEEV subroutine
        ALLOCATE(hamtot%eigval(hamtot%row), work(lwmax), rwork(3*hamtot%row-2), STAT=info)

        !CHECKPOINT: checking the allocation of the arrays
        CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it.')

        !Defining lwrok=-1 in order to compute the optimal 'lwork'
        !for the computation of the eigenvalues
        lwork=-1

        !Computing the optimal 'lwork'
        CALL ZHEEV('V', 'U', hamtot%col, hamtot%mat, hamtot%row, hamtot%eigval, work, &
             & lwork, rwork, info)

        !CHECKPOINT: checking the work of the ZHEEV subroutine
        CALL checkpoint(info/=0, 'Something wrong with the first ZHEEV. Please, check it')

        !Defining of the optimal 'lwork'
        lwork = MIN(lwmax, INT(work(1)))

        !Deallocation of the 'work' array
        DEALLOCATE(work, STAT=info)

        !CHECKPOINT: checking the deallocation
        CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please check it.')

        !Allocation of the 'work' array of the right size
        ALLOCATE(work(MAX(1, lwork)), STAT=info)

        !CHECKPOINT: checking the allocation
        CALL checkpoint(info/=0, 'Something wrong with the allocation. Please, check it.')

        !Computing of the eigenvalues
        CALL ZHEEV('N', 'U', hamtot%col, hamtot%mat, hamtot%row, hamtot%eigval, work, &
             & lwork, rwork, info)

        !CHECKPOINT: checking the work of the ZHEEV subroutine
        CALL checkpoint(info/=0, 'Something wrong with the second ZHEEV. Please, check it.')

        !Deallocation of the arrays used in the computation
        DEALLOCATE(work, rwork, hamx%mat, hamz%mat, STAT=info)

        !CHECKPOINT: checking the deallocation of the arrays
        CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it.')
     ENDIF

     !CHECKPOINT: checking the deallocation
     CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it.')

     !Storing the firs k levels
     DO jj=1, kappa
        eigenvalues(ll, jj) = hamtot%eigval(jj)
     ENDDO

     !Deallocation of the total hamiltonian and the array of eigenvalues
     DEALLOCATE(hamtot%mat, hamtot%eigval, STAT=info)

     !CHECKPOINT: checking the deallocation
     CALL checkpoint(info/=0, 'Something wrong with the deallocation. Please, check it.')


  ENDDO

  !Storing the end time of the computation
  CALL cpu_time(stopp)

  !Defining the total time of computation
  start=stopp-start

  !If statement for the printing of the right time unit
  IF (start<60.0) THEN
     PRINT*, 'Time:  ', start, 'seconds'
  ELSEIF(start>60.0 .AND. start<3600.0)THEN
     PRINT*, 'Time:  ', start/60.0, 'minuts'
  ELSE
     PRINT*, 'Time:  ', start/3600.0, 'hours'
  ENDIF

  !Writing the firs k levels in the file .gnu
  DO ii=1, kappa
     DO jj=1, binlambda
        WRITE(1, *) dl*(jj-1), eigenvalues(jj, ii)
     ENDDO
     WRITE(1,'("e")')
  ENDDO

  !Execution of the .gnu file using gnuplot
  CALL execute_command_line("gnuplot eigval.gnu -p")

  !Deallocation of all the remaining matrix
  DEALLOCATE(eigenvalues, sigmax%mat, sigmaz%mat, ide%mat, STAT=info)

  !CHECKPOINT: checking the deallocation
  CALL checkpoint(info/=0, 'Somentihg wrong with the deallcoation. Please, check it.')
  
END PROGRAM
   
