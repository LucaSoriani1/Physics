PROGRAM quantumho

   IMPLICIT NONE
   REAL*8 :: min, max, interv, dx, abs, mean, minimo, de
   REAL*8, DIMENSION(:,:), ALLOCATABLE ::  yi
   REAL*8, DIMENSION(:), ALLOCATABLE :: diag, subdiag, eigv, work,pot
   INTEGER*4, DIMENSION(:), ALLOCATABLE ::  iwork, ifail
   INTEGER*4 :: i,j, step, n, info, m
   CHARACTER :: cho
   CHARACTER(40) :: name
   DO
      PRINT*, 'What is the range?'
      PRINT*, 'Minimu:'
      READ*, min
      PRINT*, 'Maximum:'
      READ*, max
      interv=max-min
      mean=(max+min)/2
      PRINT*, 'How many small want you the discretization?'
      READ*, dx
      step=INT(interv/dx)
      step=step+1
      abs=2*0.0000000000001
      
      ALLOCATE(diag(step), subdiag(step-1), eigv(step), work(5*(step)),&
           iwork(5*(step)), ifail(step),pot(step), STAT=info) 
 
      IF (info/=0) THEN
         PRINT*, 'ALLOCATION WRONG'
         STOP
      ENDIF
 
      
      DO i=1, step
         pot(i) = ((min-mean) + (i-1)*dx)**2
      ENDDO
 
      
      DO i=1, step
         diag(i) = pot(i) + 2d0/(2*dx)**2
      ENDDO
 
      DO i=1, step-1
         subdiag(i) = -1d0/(2*dx)**2
      ENDDO
 
 
      PRINT*, 'How many eigenvalues do you want to calculate?'
      READ*, n
 
      ALLOCATE(yi(step,n))
      call dstevx('v' ,'i',step ,diag ,subdiag, 0d0,1d0, 1,n, abs, m,eigv,yi,step,&
           work,iwork,ifail,info)
   
      IF(info/=0) THEN
         PRINT*, 'The calculate fails'
         STOP
      ENDIF
      
      PRINT*, 'DO you want to store the data?'
      READ*, cho
      IF (cho=='y') THEN
         PRINT*,'In which file save the data?'
         READ*, name
         name = TRIM(name)
         OPEN(unit=1, file=name, status='unknown', action='write')
         WRITE(1, *) 'MINIMUM:  ', min
         WRITE(1, *) 'MAXIMUM:  ', max
         WRITE(1, *) 'INTERVAL:  ', interv
         WRITE(1, *) 'DISCRETIZATION:   ', dx
         WRITE(1, *) 'N. OF EIGENVALUES:   ', n
         WRITE(1, *) ''
         WRITE(1, *) 'The eigenvalues are:'
         DO i = 1, n
            WRITE(1, *) i, ':   ', eigv(i)
         ENDDO
      ENDIF
      
      DO i=1, n
         PRINT*, i,':   ', eigv(i)
      ENDDO
   
   
      CLOSE(1)
 
      PRINT*, 'Do you want use Gnuplot?'
      READ*, cho
      IF (cho=='y') THEN
         
         minimo = 0.0
         DO i=1,step
       IF(pot(i).LT.minimo) minimo=pot(i)
         ENDDO
   
         PRINT*, 'How many eigenfunction do you want to plot?'
         READ*, n
 
    de=(eigv(n)-eigv(1))/DBLE(n-1)
  
    OPEN(555,file='harmonicosci.gnu')
 
         WRITE(555,'("set nokey")')
    WRITE(555,'("set title ''Eigenfunction of harmonic oscillator in 1D''")')
    WRITE(555,'("set ylabel ''Energies''")')
    WRITE(555,'("set xlabel ''Position''")')
         WRITE(555,'("set yrange [",F12.8,":",F12.8,"]")') minimo,eigv(n)+de
         WRITE(555,'("set xrange [",I0,":",I0,"]")') INT(min), INT(max)
 
    de=1/SQRT(dx)
    WRITE(555,'("plot ''-'' w l, ",A1)') CHAR(92)
         DO i=1,n-1
              write(555,'("''-'' u 1:(",F12.8,"+",F12.8,"*($2)) w l,",A1)') eigv(i),de,char(92)
    ! WRITE(555,'("''-'' u 1:(",F12.8,"+",F12.8,"*($2)) w l,",A1)') eigv(i),de,CHAR(92)I0
    ENDDO
    WRITE(555,'("''-'' u 1:(",F12.8,"+",F12.8,"*($2)) w l")') eigv(n),de
 
    WRITE(555,'("#For plotting potential")')
    DO i=1,step
       WRITE(555,'(E20.10," ",E20.10)') min+(i-1)*dx,pot(i)
    ENDDO
    WRITE(555,'("e")')
    DO j=1,n
       WRITE(555,'("# avl(",I0,")=",E20.10)') j, eigv(j)
       DO i=1,step
          WRITE(555,'(E20.10," ",E20.10)') min+(i-1)*dx,yi(i,j)
       ENDDO
       WRITE(555,'(I0," ",E20.10)') INT(min),yi(1,j)
       WRITE(555,'("e")')
    ENDDO
    CLOSE(10)
 
    CALL execute_command_line ("gnuplot harmonicosci.gnu -p", exitstat=info)
      ENDIF
 
      DEALLOCATE(diag, subdiag, eigv, work, iwork, ifail, yi, pot)
 
 
      PRINT*, 'Do you want continue?'
      READ*, cho
      IF (cho /= 'y') THEN
         STOP
      ENDIF
      CLOSE (555)
   ENDDO
   
    
 END PROGRAM quantumho
 