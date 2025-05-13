! RSE

PROGRAM MAIN

    USE OMP_LIB

    IMPLICIT NONE

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vectorArg
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vectorSin

    INTEGER :: nList
    INTEGER, ALLOCATABLE, DIMENSION(:) :: indexList    
    
    INTEGER, PARAMETER :: N = 10240 * 10240
    
    INTEGER :: I

    DOUBLE PRECISION :: startTimer, finishTimer
    
    WRITE(*,*) 'PROGRAM START'

    ALLOCATE( vectorArg(1:N), vectorSin(1:N), indexList(1:N) )

    DO I = 1, N
        vectorArg(I) = I * (6.28 / N)
        vectorSin(I) = 0.0;
    END DO
    
    nList = 0
    indexList = -1

    startTimer  = omp_get_wtime()

    DO I = 1, N
        vectorSin(I) = SIN(vectorArg(I))
        
        IF(vectorSin(I) < 0.0) THEN
            nList = nList + 1
            indexList(nList) = I
        END IF
    END DO

    finishTimer = omp_get_wtime()
    
    WRITE(*,*) 'TIME: ', finishTimer - startTimer
    WRITE(*,*) 'LIST LENGTH: ', nList
    WRITE(*,*) 'INDEXES: ', MINVAL(indexList(1:nList)), '<> ', &
                            MAXVAL(indexList(1:nList))

     WRITE(*,*) 'PROGRAM FINISH'    

END PROGRAM MAIN
