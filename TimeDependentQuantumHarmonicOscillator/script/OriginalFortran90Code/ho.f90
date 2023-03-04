!=====================================================================!
!This program wants to calculate the eolution of the egienvectors
!of a quantum harmonic oscillator dipending on time.
!In order to solve the Schroedinger equation, I need to discretize
!the space and the time to get the right answer.
!The program consists in two part: the first is related to find the
!first eigenvalue and eigenvector of the quantum harmonic oscillator
!time indipendent, and the secondo consinst in the evolution of the
!eigenvector step by step.
!=====================================================================!




!=====================================================================!
!This subroutine compute the norm of a given complex vecotr. It takes
!in input the dimension of that vector in order to allocate it and
!as input and output the vecotr. Another input is the discretization
!in order to normalize the eigenvector. As output, the subroutine gives the
!normalized complex vector.
SUBROUTINE norm(vec, dim, dx)
  
  IMPLICIT NONE
  
  !Dimension of the input vector to allocate it
  INTEGER*8, INTENT(IN) :: dim

  !The input/output allocated complex vector
  COMPLEX*16, INTENT(INOUT):: vec(dim)

  !Discretization
  REAL*8, INTENT(IN) :: dx

  !The real the will be the norm of the vector
  REAL*8 :: norma

  !Index for the DO cycle
  INTEGER*8 :: i

  !Inizialized of the norm to zero
  norma = 0d0

  !Computation of the norm
  Do i=1, SIZE(vec)
     norma = norma + (REAL(vec(i))**2 + AIMAG(vec(i))**2)*dx
  ENDDO
  norma = SQRT(norma)
  vec = vec/norma
  
  !Checking if the norm was right
  norma=0d0
  Do i=1, SIZE(vec)
     norma = norma + (REAL(vec(i))**2 + AIMAG(vec(i))**2)*dx
  ENDDO
  CALL checkpoint(ABS(norma-1.0) > 0.1, 'SOMETHING WAS WRONG WITH THE NORM')
  !Normalization of the vector
END SUBROUTINE norm
!=====================================================================!



!=====================================================================!
!This subroutine checks if a condition is true or false. As input it
!take the condition and a message, that is the error message. In the
!main program, the condition is represented by the wrong condition: if
!this wrong condition is true, the program will print the error message
!given as input and will print the position of the check in order to
!simplify the search of the problem.
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
!=====================================================================!
!--------------------------MAIN PROGRAM-------------------------------!
!=====================================================================!
!=====================================================================!
PROGRAM quantumhotime

  IMPLICIT NONE
  
  !This file is necessary to use the fftw3 and it contains a series of
  !parameters that fftw3 will use
  INCLUDE 'fftw3.f'

  !******************************************************************!
  !Index for the do cycles
  INTEGER*4 :: i, j, index

  !******************************************************************!
  !Spatial variables
  !Default parameters for the spatial range
  REAL*8 :: min=-15.0, max=15.0, dx=0.01 !minimum and maximum of the 
                                          !the range and the
                                          !discretization of the space  
  INTEGER*8 :: intspace=30, binspace  !total range and bin of the space

  !******************************************************************!  
  !First eigenvalue variables: array to create the Hamiltonian
  REAL*8, ALLOCATABLE :: diag(:), subdiag(:) !diagonal and subdiagonal

  !******************************************************************!
  !Potential variables: array used to store the potential of the HO
  !the kinetic part
  REAL*8, ALLOCATABLE :: potx(:), potk(:)

  !******************************************************************!
  !DSTEVX variables: the first array and the first matrix are used to
  !store respectively the first eigenvalue and the eigenvector
  REAL*8, ALLOCATABLE :: eigenvalues(:), eigenvectors(:,:), work(:)
  REAL*8 :: abstol=2*0.0000000000001 
  INTEGER*4 :: m, info
  INTEGER*4, ALLOCATABLE :: iwork(:), ifail(:)

  !******************************************************************!
  !Time variables
  INTEGER*4 :: T=100, bintime !Total time and bin for the time
  REAL*8 :: dt = 0.01, start, end !Discretization and variables used
                                  !to check the computation time
  
  REAL*8, ALLOCATABLE :: time(:) !Array used to store each time-step
  COMPLEX*8 :: ii = CMPLX(0.0, 1.0) !COMPLEX UNIT!

  !******************************************************************!
  !Time-evolution variables
  !The first is the array to store the eigenvector, the second is
  !used to store a copy of the first. The third is the vector where
  !the program stores the discretized momenta of the particle. The
  !fourth is used to stored the Fourier transform of the 'temp' array
  COMPLEX*16, ALLOCATABLE :: psi(:), temp(:), k_i(:), temp1(:)

  !Variable used to discretize the momenta
  REAL*8 :: dk, k, pi=3.14159265358979, sum

  !******************************************************************!  
  !FFTW variables: is store the FORWARD subroutine and the BACKWARD
  !subroutine in these two integers
  INTEGER*8 :: for, back

  !******************************************************************!
  !Mean of the position
  REAL*8, ALLOCATABLE :: means(:)
  REAL*8 :: somma
  INTEGER*8 :: top(1)

  !******************************************************************!
  !Choise variable
  CHARACTER(6) :: cho
  INTEGER :: scelta

  !******************************************************************!
  !******************************************************************!
  !-------------------------INITIALIZATION---------------------------!
  !******************************************************************!
  !******************************************************************!

  
  !******************************************************************!
  !In this part the user can see the default parameters and, if he
  !wants, can change these parameters, by writing "change"
  PRINT*, '============================================================================'
  PRINT*, 'TIME EVOLUTION OF QUANTUM HARMONIC OSCILLATOR WITH TIME-DEPENDET HAMILTONIAN'
  PRINT*, '============================================================================'
  PRINT*,'The parameters are:'
  WRITE(*,'("---> L=",I0, " (Min=", I0, ", Max=", I0,")")') intspace, INT(min), INT(max)
  WRITE(*, '("---> dx = ", F5.3 )') dx
  WRITE(*, '("---> T = ", I0 )') T
  WRITE(*, '("---> dt = ", F4.2 )') dt
  PRINT*, 'Type "change" if you want to change them or type something else to continue'
  READ*, cho
  IF (cho=='change') THEN
     !This DO cycles if infinite once the user type 0. This is made
     !in order to allow the user to change all this parameters as he
     !whishes
     DO
        PRINT*,''
        PRINT*, '============================================================================'
        PRINT*, 'What you want to change? (Press 0 if you want to exit)'
        WRITE(*,'(" 1) L=",I0, " (Min=", I0, ", Max=", I0,")")') intspace, INT(min), INT(max)
        WRITE(*, '(" 2) dx = ", F7.5 )') dx
        WRITE(*, '(" 3) T = ", I0 )') T
        WRITE(*, '(" 4) dt = ", F6.4 )') dt
        PRINT*,'0) EXIT'
        READ*, scelta
        IF (scelta==1) THEN           
           PRINT*,'Set the minimum'
           READ*, min
           PRINT*, 'set the maximum'
           READ*, max
           !POSITION #0.1: checking if the user inserts the same values
           CALL checkpoint( max==min, 'YOU MUST INSERT DIFFERENT VALUES')
           intspace = max-min
           PRINT*, '============================================================================'
           PRINT*, ''
           WRITE(*, '("The range is L=", I0)') intspace
           PRINT*, ''
        ELSE IF (scelta == 2) THEN
           PRINT*, 'Set the discretization in space'
           READ*, dx
           !POSITION #0.2: checking if the discretization is right
           CALL checkpoint( dx <= 0 .OR. dx >=1, 'YOU MUST INSERT A dx BETWEEN 0 AND 1')
        ELSE IF (scelta == 3) THEN
           PRINT*, 'Set the maximum Time'
           READ*, T
           !POSITION #0.3: checking if the time is right
           CALL checkpoint (T<=1, 'YOU MUST INSERT A TIME BIGGER THEN 1')
        ELSE IF (scelta == 4) THEN
           PRINT*, 'Set the discretizaztion in time'
           READ*, dt
           !POSITION #0.4: checking if the discretization is right
           CALL checkpoint( dt <= 0 .OR. dt >=1, 'YOU MUST INSERT A dt BETWEEN 0 AND 1')
        ELSE IF (scelta == 0) THEN
           EXIT
        ELSE
           PRINT*, 'You make no choice'
           PRINT*, ''
           
        ENDIF
     ENDDO
     !This is 'secret' key-word in order to exit from this IF statement
  ELSE IF (cho=='exittt') THEN
     STOP
  ENDIF

  !Defining the bins for space and time with the default/user parameters
  bintime = T/dt
  binspace = INT(intspace/dx)+1   

  !Definition of the discretization of the momenta
   dk = pi/intspace
   !POSITION #1: checking if dk is bigger than 0
   CALL checkpoint( dk <= 0, 'ERROR IN THE COMPUTATION OF dk. PLEASE RESTART THE PROGRAM')
  
  !Allocatation of the vectors
  ALLOCATE(potx(binspace), diag(binspace), subdiag(binspace-1), eigenvalues(binspace), &
       eigenvectors(binspace,1), work(5*binspace), iwork(5*binspace), ifail(binspace), &
       psi(binspace), potk(binspace), temp(binspace),k_i(binspace), means(bintime+1), &
       time(bintime),temp1(binspace), STAT=info)

  !POSITION #2: checking the allocation
  CALL checkpoint(info/=0, 'ALLOCATION WAS WRONG, PLEASE RESTART')
  
  !******************************************************************!
  !******************************************************************!
  !--------------------BEGINNING OF THE CALCULATION------------------!
  !******************************************************************!
  !******************************************************************!


  
  !******************************************************************!
  !-------------------FIRST PART: EIGENPROBLEM-----------------------!
  !******************************************************************!

  !Storing the time step
  DO i=1, bintime+1
     time(i) = (i-1)*dt/T
  ENDDO
  
  
  !Defining the potential of the time-indipendet HO Hamiltonian
  !N.B. the 1/2 factor is not necessary because the Hamiltonian
  !is indipendent by scaling factor
  DO i=1, binspace
     potx(i) = (min + (i-1)*dx)**2
  ENDDO

  !******************************************************************!
  !The discretized Hamiltonian is formed by 2 on the main diagonal and
  ! by -1 on the two sub-diagonal (tridiagonal matrix). I need to
  !divide these values by 2 times the spatial discretization dx because
  !we obtained the Hamiltonian by computing the central derivative.
  !******************************************************************!
  
  !Defining the diagonal elements
  DO i=1, binspace
     diag(i) = potx(i) + (2d0/(2*dx)**2)
  ENDDO

  !Defining the sub-diagonal elements
  DO i=1, binspace-1
     subdiag(i) = -1d0/(2*dx)**2
  ENDDO

  !Calculation of the first eigenvalue and eigenvector
  CALL dstevx('v', 'i', binspace, diag, subdiag, 0d0, 1d0, 1, 1, abstol, m, eigenvalues, &
       eigenvectors, binspace, work, iwork, ifail, info)
  !POSITION #3: checking the work of DSTEVX
  CALL checkpoint( info/=0, 'SOMETHING WAS WRONG WITH THE COMPUTATION OF EIGENVALUES')
  
  !POSITION #4: checking the right values of the eigenvalue
  CALL checkpoint(ABS(eigenvalues(1)-0.5)>0.1, 'SOMETHING WAS WRONG WITH THE EIGENVALUE, TRY TO SET DIFFERENT PARAMETERS')
  print*, ''
  PRINT*, '============================================================================'
  PRINT*, 'The first eigenvalue is:  ', eigenvalues(1)
  PRINT*, 'Type anything to continue'
  READ*, cho
  !'Secret' key-word to exit from the program
  IF (cho == 'exittt') THEN
     STOP
  ENDIF
  
  
  !******************************************************************!
  !------------------SECOND PART: TIME EVOLUTION---------------------!
  !******************************************************************!

  !Defining psi=first_eigenvector: psi is a complex where, for this
  !moment is just formed by the real part. I divide by the square root
  !of the discretization in order to normalize the eigenvector
  DO i=1, binspace
     psi(i) = CMPLX(eigenvectors(i,1), 0d0)
  ENDDO

  !Normalizing the vector
  CALL norm(psi, binspace, dx)

  
  !******************************************************************!
  !The outputs of the fftw are order in the following way: the first
  !part of the array refers to that momenta with positive frequencies
  !and the second parts refers to that momenta with negative
  !fequencies. Thus I'have ordered the momenta following this sets
  !******************************************************************!
  
  !Storing the momenta for the Fourier transform in according with
  !the output of the fftw3
  DO i=2, (SIZE(k_i)+1)/2
     k_i(i)=(i-1)*dk
     k_i(SIZE(k_i)+2-i) =(i-1)*dk
  ENDDO
  k_i((SIZE(k_i)+2)/2) = pi/dx
  
  !******************************************************************!
  !OPEN AND PRINTING stuff on the .gnu file
  OPEN(unit=10, file='evol1.gnu')
  WRITE(10, '("set nokey")')
  WRITE(10, '("set grid")')
  WRITE(10, '("set title ''H.O. with time-dependet Hamiltonian evolution (L=",I0,", dx=",F6.4,", T=",I0,", dt=",F5.3,")''")') &
       intspace, dx, T, dt
  WRITE(10, '("set ylabel ''Square modulus of Psi''")')
  WRITE(10, '("set xlabel ''Position''")')
  WRITE(10, '("set xrange [", I0,":",I0,"]")') INT(min/4), INT(max/4)+1
  WRITE(10, '("set yrange [0:1]")')
  WRITE(10,'("plot ''-'' u 1:2 w l title ''ground state'',",A1)')CHAr(92)
  DO i=1,bintime-1
     WRITE(10,'(" ''-'' u 1:2 w l title ''q_t=",ES9.3,"'',",A1)')DBLE(i/T)*dt,CHAR(92)
     WRITE(10,'(" ''-'' u 1:2 w l title ''pot_t=",ES9.3,"'',",A1)')DBLE(i/T)*dt,CHAR(92)
  ENDDO
  WRITE(10,'(" ''-'' u 1:2 w l title ''q_t=",ES9.3,"'',",A1)')DBLE((i+1)/T)*dt,CHAR(92)
  WRITE(10,'(" ''-'' u 1:2 w l title ''pot_=",ES9.3,"''")')DBLE((i+1)/T)*dt
  WRITE(10,'("# primo autovettore ")')
  !******************************************************************!


  !Computation of the mean position for the ground_state
  somma=0d0
  DO i =1, binspace
     somma = somma + dx*dx*(REALPART(psi(i))**2 + AIMAG(psi(i))**2)*(i-1)
  ENDDO
  means(1) = somma
  
  !Writing in the file the ground_state
  DO i=1, binspace
     WRITE(10,*) min+(i-1)*dx, REAL(psi(i))**2+AIMAG(psi(i))**2
  ENDDO
  WRITE(10,'("e")')

  
  !******************************************************************!
  !******************************************************************!
  !-------------------------TIME EVOLUTION---------------------------!
  !******************************************************************!
  !******************************************************************!
  
  !Index used to calculate the percentage of computation
  index = 0

  !Variable used to calculate the time of computation
  CALL cpu_time(start)
  DO j=1, bintime
     !IF statement for the percentage
     IF ( INT(DBLE(j)/DBLE(bintime)*100d0) == index) THEN
        WRITE(*, '("---> ",I0, "%")')index
        index = index +1
     ENDIF
     !Copying the psi in a dummy array
     temp=psi
     
     !Defining the potenzial depending on time
     potx=0d0
     DO i=1, binspace
        potx(i) = ((min + (i-1)*dx - j*dt/T)**2)/2
     ENDDO

     !Temporal evolution of the eigenvector
     DO i=1, binspace
        temp(i) = EXP(-ii* potx(i) * dt/2)* temp(i)
     ENDDO
     
    !Fourier transform of the temporal evoluted eigenvector
     CALL dfftw_plan_dft_1d(for, binspace, temp, temp1, FFTW_FORWARD,FFTW_ESTIMATE)
     CALL dfftw_execute_dft(for, temp, temp1)
     CALL dfftw_destroy_plan(for)
     
     !Normalizing the output eigenvector
     CALL norm(temp1, binspace, dk)

     !Defining the kinetic part of the Hamiltonian
     DO i=1, binspace
        potk(i) = (k_i(i)**2)/2
     ENDDO
     
     !Temporal evolution of the tranformed eigenstate
     Do i=1, binspace
        temp1(i) = EXP(-ii*potk(i)*dt)*temp1(i)
     ENDDO

     !Anti-Fourier tranform
     CALL dfftw_plan_dft_1d(back, binspace, temp1, temp, FFTW_BACKWARD,FFTW_ESTIMATE)
     CALL dfftw_execute_dft(back, temp1, temp)
     CALL dfftw_destroy_plan(back)

     !Normalizing the output eigenvector
     CALL norm(temp, binspace, dx)

     !Last temporal evolution
     DO i=1, binspace
        temp(i) = EXP(-ii*potx(i)*dt/2)*temp(i)
     ENDDO

     !Re-normalizing the eigenvector
     CALL norm(temp, binspace, dx)
     
     !Assing to the initial psi the evolved eigenvector
     psi=temp

     !Writing the square module in the file
     DO i=1, binspace
        WRITE(10, *) min+(i-1)*dx, REAL(psi(i))**2 + AIMAG(psi(i))**2
     ENDDO   
     WRITE(10, '("e")')

     !Writing the evolved potential in the file
     DO i=1,binspace
        WRITE(10,*) min+(i-1)*dx, potx(i)
     ENDDO
     WRITE(10, '("e")')

     !Storing the mean position
     somma = 0d0
     DO i =1, binspace
        somma = somma + dx*dx*(REAL(psi(i))**2 + AIMAG(psi(i))**2)*(i-1)
     ENDDO
     means(j+1)=somma
  ENDDO
  !Calculation of the time spent in minutes
  CALL cpu_time(end)
  start=end-start
  PRINT*, 'I spent ', start/60, 'minutes to do this!'

  !Execution of the output 
  CALL execute_command_line("gnuplot evol1.gnu -p")

  !Closing the fil
  CLOSE(10)

  
  !******************************************************************!
  !Opening and creating a new .gnu file in order to plot the mean position
  OPEN(unit=3, file='timeevol.gnu')
  WRITE(3, '("set grid")')
  WRITE(3, '("f(x) = x")')
  WRITE(3, '("set title ''Evolution of the position (L=",I0,", dx=",F6.4,", T=",I0,", dt=",F5.3,")''")')&
       intspace, dx, T, dt
  WRITE(3, '("set ylabel ''Mean position''")')
  WRITE(3, '("set xlabel ''TIme step''")')
  WRITE(3,'("plot ''-'' u 1:2 w l title ''Position at each time-step'', f(x) title ''Minimum of the potential''")')
  DO i=1, bintime+1
     WRITE(3, *) time(i), means(i)+min
  ENDDO
  CLOSE(3)
  !******************************************************************!

  !Execution of the .gnu file
  CALL execute_command_line("gnuplot timeevol.gnu -p")
  
  DEALLOCATE(potx, diag, subdiag, eigenvalues, eigenvectors, work,&
       iwork, ifail, psi, potk, temp, k_i, means, time, STAT=info)
  CALL checkpoint(info/=0, 'SOMETHING WAS WRONG IN THE DEALLOCATION')

END PROGRAM quantumhotime
