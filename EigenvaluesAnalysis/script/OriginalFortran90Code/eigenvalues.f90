!gfortran -L/usr/lib/x86_64-linux-gnu/lapack -o exe5 exe5_1.f90 -llapack

MODULE matrix
  TYPE dmatrix
     INTEGER*4 :: ddim
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: dmatrice
  END TYPE dmatrix

  TYPE dmatrix2
     INTEGER*4 :: ddim2
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: dmatrice2
  END TYPE dmatrix2

  TYPE hermitian
     INTEGER*4 :: ddimh
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: dmatriceh
  END TYPE hermitian
  
  INTERFACE OPERATOR (.create.)
     MODULE PROCEDURE creatematrix
  END INTERFACE OPERATOR(.create.)

  INTERFACE OPERATOR (.adjoint.)
     MODULE PROCEDURE adjointmatrix
  END INTERFACE OPERATOR (.adjoint.)

  INTERFACE OPERATOR (.hermi.)
     MODULE PROCEDURE hermitianmatrix
  END INTERFACE OPERATOR (.hermi.)

  INTERFACE OPERATOR (.U.)
     MODULE PROCEDURE triangularU
  END INTERFACE OPERATOR (.U.)
  
  INTERFACE OPERATOR (.L.)
     MODULE PROCEDURE triangularL
  END INTERFACE OPERATOR (.L.)
  
  
  
  
  
  
CONTAINS
  TYPE(dmatrix) FUNCTION creatematrix(dim)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: dim
    INTEGER*4 :: i,j
    REAL*16, DIMENSION(:,:), ALLOCATABLE :: re, im
    REAL*16 :: a,b

    ALLOCATE(re(dim, dim), im(dim, dim))
    ALLOCATE(creatematrix%dmatrice(dim, dim))

    creatematrix%ddim = dim

    CALL random_number(re)
    CALL random_number(im)

    DO i=1, dim
       DO j=1, dim
          IF (i==j) THEN
             CALL random_number(a)
             a = 2*a-1
             creatematrix%dmatrice(i,j) = CMPLX(a, 0.0)
          ELSE
             CALL random_number(a)
             CALL random_number(b)
             a=2*a-1
             b=2*b-1
             creatematrix%dmatrice(i,j) = CMPLX(a,b)
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(re, im)
  END FUNCTION creatematrix

  TYPE(dmatrix2) FUNCTION adjointmatrix(mat2)
    
    IMPLICIT NONE
    INTEGER*4 :: i,j
    TYPE(dmatrix), INTENT(IN) :: mat2

    adjointmatrix%ddim2 = mat2%ddim
    
    ALLOCATE(adjointmatrix%dmatrice2(adjointmatrix%ddim2, adjointmatrix%ddim2))
    DO i=1, adjointmatrix%ddim2
       DO j=1, adjointmatrix%ddim2
          adjointmatrix%dmatrice2(i,j) = CMPLX(REAL(mat2%dmatrice(j,i)), -AIMAG(mat2%dmatrice(j,i)))
       ENDDO
    ENDDO

  END FUNCTION adjointmatrix

  TYPE(hermitian) FUNCTION hermitianmatrix(mat1, mat2)

    IMPLICIT NONE
    TYPE(dmatrix), INTENT(IN) :: mat1
    TYPE(dmatrix2), INTENT(IN) :: mat2

    hermitianmatrix%ddimh = mat1%ddim
    

    ALLOCATE(hermitianmatrix%dmatriceh(hermitianmatrix%ddimh, hermitianmatrix%ddimh))

    hermitianmatrix%dmatriceh = 0.5*(mat1%dmatrice + mat2%dmatrice2)
    

  END FUNCTION hermitianmatrix

  TYPE(hermitian) FUNCTION triangularU(mat2)
    
    IMPLICIT NONE
    INTEGER*4 :: i,j
    TYPE(hermitian), INTENT(IN) :: mat2

    triangularU%ddimh = mat2%ddimh

    ALLOCATE(triangularU%dmatriceh(mat2%ddimh, mat2%ddimh))
    
    DO i=1, mat2%ddimh
       DO j=1, mat2%ddimh
          IF (i<=j) THEN
             triangularU%dmatriceh(i,j) = mat2%dmatriceh(i,j)
          ELSE
             triangularU%dmatriceh(i,j)=CMPLX(0.0,0.0)
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION triangularU

  
  TYPE(hermitian) FUNCTION triangularL(mat2)

    IMPLICIT NONE
    INTEGER*4 :: i,j
    TYPE(hermitian), INTENT(IN) :: mat2

    triangularL%ddimh = mat2%ddimh

    
    ALLOCATE(triangularL%dmatriceh(mat2%ddimh, mat2%ddimh))
    
    DO i=1, mat2%ddimh
       DO j=1, mat2%ddimh
          IF (i>=j) THEN
             triangularL%dmatriceh(i,j) = mat2%dmatriceh(i,j)
          ELSE
             triangularL%dmatriceh(i,j)=CMPLX(0.0,0.0)
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION triangularL
  
    
    
  
END MODULE matrix


PROGRAM principale

  USE matrix
  IMPLICIT NONE
  INTEGER*4 :: dim, l, max, i,j,info, lwmax=10000, lwork
  TYPE(dmatrix) :: matrice1
  TYPE(dmatrix2) :: matrice2
  TYPE(hermitian) :: matriceh, triang1, triang
  REAL*8, DIMENSION(:), ALLOCATABLE :: w,rwork, medie, spa, w2, medie2, spa2, w3
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: diagon
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: work
  CHARACTER :: tri,cho, pri
  CHARACTER*40 :: name, name2
  REAL*8 :: mean, time1, time2, end, mini, maxi, step, x

  PRINT*, 'What is the values of the dimension N?'
  READ*, dim
  PRINT*, 'How many cycles do you want to do?'
  READ*, max
  PRINT*, 'In which file do you want store the data for the Hermitian matrix?'
  READ*, name
  PRINT*, 'In which file do you want store the data for the Hermitian matrix?'
  READ*, name2
  name=TRIM(name)
  name2= TRIM(name2)
  CALL cpu_time(time2)
  OPEN(1, file=name, status='unknown', position='append', action='write')
  OPEN(2, file=name2, status='unknown', position='append', action='write')
  DO l=1, max
     CALL cpu_time(time1)
     PRINT*, ''
     PRINT*, l
     matrice1 = .create.dim
     matrice2 = .adjoint.matrice1
     matriceh = matrice1.hermi.matrice2
     tri='U'
     cho='N'  

     triang = .U.matriceh
     triang1 = triang
     ALLOCATE(w(dim), rwork(3*dim-2),work(lwmax))
     lwork=-1
     CALL ZHEEV(cho, tri, dim, triang%dmatriceh, dim, w,work,lwork,rwork,info)
     
     lwork=MIN(lwmax, INT(work(1)))
     CALL ZHEEV(cho, tri, dim, triang%dmatriceh, dim, w,work,lwork,rwork,info)  

     ALLOCATE(w2(SIZE(w)-1))

     DO i =1, SIZE(w2)
        w2(i)=w(i)
     ENDDO

     ALLOCATE(medie(dim-2),spa(dim-2))
     mean=0.0
     DO i=1, SIZE(spa)
        medie(i) = w2(i+1) - w2(i)
        mean = mean + medie(i)
     ENDDO
     mean=mean/SIZE(medie)
     DO i=1, dim-2
        spa(i) = medie(i)/mean
        WRITE(1,*) spa(i)
     ENDDO

     
     DEALLOCATE(w, rwork, work)
     ALLOCATE (w(dim), rwork(3*dim-2), work(lwmax))
     lwork=-1
     
     ALLOCATE(diagon(dim, dim))
     DO i=1, dim
        DO j=1, dim
           IF(i==j)THEN
              CALL random_number(x)
              x=2*x-1
              diagon(i,j) = CMPLX(x, 0.0)
           ELSE
              diagon(i,j)=CMPLX(0.0,0.0)
           ENDIF
        ENDDO
     ENDDO

     CALL ZHEEV(cho, tri, dim, diagon, dim, w, work, lwork, rwork, info)
     lwork = MIN(lwmax, INT(work(1)))
     CALL ZHEEV(cho, tri, dim, diagon, dim, w, work, lwork, rwork, info)


     ALLOCATE(medie2(dim-1), spa2(dim-1))
     mean=0.0
     DO i=1, SIZE(spa2)
        medie2(i) = w(i+1)-w(i)
        mean= mean+medie2(i)
     ENDDO
     mean=mean/SIZE(medie2)

     DO i=1, SIZE(medie2)
        spa2(i) = medie2(i)/mean
        WRITE(2,*) spa2(i)
     ENDDO
    

     CALL cpu_time(end)
     time1=end-time1
     PRINT*, 'TIME: ', time1, 's'
     
     DEALLOCATE(w, rwork, work, diagon, medie,w2,spa, medie2, spa2)
  ENDDO
  CALL cpu_time(end)
  time2=end-time2
  PRINT*, 'TOTAL TIME: ', time2/60, ' minuts'
  CLOSE(1)           
  CLOSE(2)
END PROGRAM principale
     
