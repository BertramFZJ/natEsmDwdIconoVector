PROGRAM MAIN

    USE OMP_LIB

    IMPLICIT NONE

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: matrixARG
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: matrixSIN2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: matrixCOS2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: matrixRES
    
    INTEGER, PARAMETER :: N = 10240
    
    INTEGER :: I, J, ID

    DOUBLE PRECISION :: startTimer, finishTimer
    
    WRITE(*,*) 'PROGRAM START'

    ALLOCATE( matrixARG(1:N,1:N), matrixSIN2(1:N,1:N), matrixCOS2(1:N,1:N), matrixRES(1:N,1:N) )

    startTimer = omp_get_wtime()

    DO J = 1, N
        DO I = 1, N
            ID = I + (J - 1) * N - 1
            matrixARG(I,J) = ID * 0.0123
        END DO
    ENDDO

    !$OMP PARALLEL DO PRIVATE(I) NUM_THREADS(4)
    DO J = 1, N
        !$OMP SIMD
        DO I = 1, N
            matrixSIN2(I,J) = SIN(matrixARG(I,J))
            matrixSIN2(I,J) = matrixSIN2(I,J) * matrixSIN2(I,J)
            matrixCOS2(I,J) = COS(matrixARG(I,J))
            matrixCOS2(I,J) = matrixCOS2(I,J) * matrixCOS2(I,J)
        END DO
    ENDDO
    !$OMP END PARALLEL DO

    DO J = 1, N
        !$OMP SIMD
        DO I = 1, N
            CALL subSUM(matrixSIN2(I,J), matrixCOS2(I,J), matrixRES(I,J))
        END DO
    ENDDO

    CALL subSUM(matrixSIN2(:,:), matrixCOS2(:,:), matrixRES(:,:))

    finishTimer = omp_get_wtime()
    WRITE(*,*) "TIME ", finishTimer - startTimer

    WRITE(*,*) 'PROGRAM FINISH'

    CONTAINS

    PURE ELEMENTAL SUBROUTINE subSUM(a, b, c)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(in)    :: a
        DOUBLE PRECISION, INTENT(in)    :: b
        DOUBLE PRECISION, INTENT(inout) :: c

        c = a + b    

    END SUBROUTINE subSUM

END PROGRAM MAIN