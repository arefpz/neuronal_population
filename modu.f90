!! This Program is the module used in the simulation of following paper

!! Compiled on 21 Feb. 2021 by gfortran 10.2.0 on Ubuntu 20.10 as following lines

!! gfortran-m64-O3 hh_mod.f95 main_hh.f95-o a1.out
!! echo $deeout $dI | ./a1.out

!! Written by Aref Pariz and Farhad Daei. Email: pariz.aref@gmail.com

MODULE modu
    USE VARS
    IMPLICIT NONE
    
CONTAINS
    SUBROUTINE initialize_seed(i)
        INTEGER:: i
        INTEGER, DIMENSION(8):: dtt
        INTEGER, DIMENSION(:), ALLOCATABLE:: seed

        CALL RANDOM_SEED(SIZE = i)
        ALLOCATE(seed(i))
        CALL DATE_AND_TIME(VALUES = dtt)  ! dtt = 1
        CALL RANDOM_SEED(GET = seed)
        seed(i) = dtt(8)
        seed(1) = dtt(8) * dtt(7) * dtt(6)
        CALL RANDOM_SEED(PUT = seed)
    END SUBROUTINE initialize_seed
    
    FUNCTION funv(vv, vr, synI, extI, taum)
        IMPLICIT NONE
        REAL(KIND = DP), INTENT(IN):: vv, synI, extI, vr
        REAL(KIND = DP):: funv, taum
        ! Euler-Maruyama Method
        funv = (vr+extI+synI  - vv) / taum
    END FUNCTION funv

    FUNCTION intrnd(a, b, c) result(intrnd_r)
        IMPLICIt none
        LOGICAL :: similar
        integer(KIND = 4) :: a, b, c
        INTEGER :: nc
        real(KIND = 4):: x
        INTEGER(KIND=4) :: y
        integer(KIND = 4), dimension(:):: intrnd_r(c)
        integer(KIND = 4) :: in 
        intrnd_r = 0
        nc = 0
        DO WHILE (nc .LT. c)
            call random_number(x)
            y = a + int((b-a) * x)
            similar = .FALSE.
            DO in = 1, c    
                IF (y .EQ. intrnd_r(in)) THEN
                    similar = .TRUE.
                    cycle
                END IF
            END DO
            IF (similar .EQV. .FALSE.) THEN
                nc = nc + 1
                    intrnd_r(nc) = y
            END IF
        END DO
    END function intrnd

    FUNCTION creatA_probability(N, Ne, Ni, Pee1, Pei1, Pie1, Pii1)
        IMPLICIT NONE
        ! INTEGER:: nPee, nPei, nPie, nPii
        !        REAL(KIND = DP), ALLOCATABLE, DIMENSION(:):: xpee, xpei, xpie, xpii
        ! INTEGER, ALLOCATABLE, DIMENSION(:):: ypee, ypei, ypie, ypii
        INTEGER, INTENT(IN):: Ne, Ni, N
        REAL(KIND = DP), INTENT(IN):: Pee1, Pei1, Pie1, Pii1
        INTEGER:: i1, j1, en1, st1, st2, en2
        INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:):: creatA_probability
        REAL :: rand_number
        ALLOCATE(creatA_probability(N, N))
        creatA_probability(1: N, 1: N) = 0
        !! E -> E
        st1 = 1
        en1 = Ne
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1 ! Rows
            DO j1 = st2, en2  ! Columns
                CALL random_number(rand_number)
                IF (rand_number .LE. Pee1) THEN
                    creatA_probability(i1, j1) = 1
                END IF
            END DO
        END DO
        !! I -> E
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2  ! Columns
                CALL random_number(rand_number)
                IF (rand_number .LE. Pie1) THEN
                    creatA_probability(i1, j1) = 1
                END IF
            END DO
        END DO
        !! E -> I
        st1 = 1
        en1 = Ne
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2
                CALL random_number(rand_number)
                IF (rand_number .LE. Pei1) THEN
                    creatA_probability(i1, j1) = 1
                END IF
            END DO
        END DO
        !! I -> I
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1
            DO j1 = st2, en2
                CALL random_number(rand_number)
                IF (rand_number .LE. Pii1) THEN
                    creatA_probability(i1, j1) = 1
                END IF
            END DO
        END DO
    END FUNCTION creatA_probability

    FUNCTION creatA_fixed_number(N, Ne, Ni, Pee1, Pei1, Pie1, Pii1)
        IMPLICIT NONE
        INTEGER:: nPee, nPei, nPie, nPii
        INTEGER(KIND = DP), ALLOCATABLE, DIMENSION(:):: ypee, ypei, ypie, ypii
        INTEGER, INTENT(IN):: Ne, Ni, N
        REAL(KIND = DP), INTENT(IN):: Pee1, Pei1, Pie1, Pii1
        INTEGER:: i1, j1, en1, st1, st2, en2
        INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:):: creatA_fixed_number
        ALLOCATE(creatA_fixed_number(N, N))
        creatA_fixed_number(1: N, 1: N) = 0
        nPee = INT(Ne*Pee1)
        nPei = INT(Ne*Pei1)
        nPie = INT(Ni*Pie1)
        nPii = INT(Ni*Pii1)
        ALLoCATE(ypee(npee),ypei(nPei),ypie(nPie),ypii(nPii))
        !! E -> E
        st1 = 1
        en1 = Ne
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1  ! Rows ! Post
            ypee = intrnd(st2, en2, nPee)
            DO j1 = 1, SIZE(ypee)  ! Columns ! Pre
                creatA_fixed_number( ypee(j1), i1) = 1
            END DO
        END DO
        !! E -> I
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1  ! Rows ! Post
            ypei = intrnd(st2, en2, nPei)
            DO j1 = 1, SIZE(ypei)  ! Columns !Pre
                creatA_fixed_number(ypei(j1), i1) = 1
            END DO
        END DO
        !! I -> E
        st1 = 1
        en1 = Ne
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1  ! Rows ! Post
            ypie = intrnd(st2, en2, nPie)
            DO j1 = 1, SIZE(ypie) ! Pre
                creatA_fixed_number(ypie(j1), i1) = 1
            END DO
        END DO
        !! I -> I
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1
            ypii = intrnd(st2, en2, nPii)
            DO j1 = 1, SIZE(ypii)
                creatA_fixed_number(ypii(j1), i1) = 1
            END DO
        END DO
    END FUNCTION creatA_fixed_number

    FUNCTION create_delay(dee, dei, die, dii, ddee, ddei, ddie, ddii, neurons, Ne, Ni, dt, A)
        implicit none
        INTEGER, INTENT(IN):: Ne, Ni, neurons
        INTEGER(KIND = 1), INTENT(IN), DIMENSION(:,:) :: A
        INTEGER :: st1, st2, en1, en2, j1, i1
        REAL(KIND = DP), INTENT(IN):: dee, dei, die, dii
        REAL(KIND = DP) :: dt
        REAL(KIND = DP), INTENT(IN):: ddee, ddei, ddie, ddii
        !INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:):: A
        INTEGER(KIND = 4), DIMENSION(:,:):: create_delay( neurons,  neurons)
        REAL(KIND = DP) :: daux

        !! E -> E
        st1 = 1
        en1 = Ne
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2  ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(daux)
                    create_delay( i1, j1) = INT((dee + ddee*(-1 + 2*daux))/dt)
                END IF
            END DO
        END DO
        !! E -> I
        st1 =  1
        en1 =  Ne
        st2 = Ne+1
        en2 = Ne+Ni
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2  ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(daux)
                    create_delay( i1, j1) = INT((dei + ddei*(-1 + 2*daux))/dt)
                END IF
            END DO
        END DO
        !! I -> E
        st1 = Ne+1
        en1 = Ne+Ni
        st2 =  1
        en2 =  Ne
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2 ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(daux)
                    create_delay( i1, j1) = INT((die + ddie*(-1 + 2*daux))/dt)
                END IF
            END DO
        END DO
        !! I -> I
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1
            DO j1 = st2, en2 ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(daux)
                    create_delay( i1, j1) = INT((dii + ddii*(-1 + 2*daux))/dt)
                END IF
            END DO
        END DO


    END FUNCTION create_delay

    FUNCTION synapses(gee1, gei1, gie1, gii1, dgee1, dgei1, dgie1, dgii1, neurons, Ne, Ni, A)
        IMPLICIT NONE
        INTEGER :: st1, st2, en1, en2, j1, i1
        INTEGER, INTENT(IN):: Ne, Ni, neurons
        REAL(KIND = DP), INTENT(IN):: gee1, gei1, gie1, gii1
        REAL(KIND = DP), INTENT(IN):: dgee1, dgei1, dgie1, dgii1
        INTEGER(KIND = 1), INTENT(IN), DIMENSION(:,:) :: A
        REAL(KIND = DP), ALLOCATABLE, DIMENSION (:,:):: synapses
        REAL :: gaux
        ALLOCATE(synapses( neurons,  neurons))
        synapses = 0
        !! E -> E
        st1 = 1
        en1 = Ne
        st2 = 1
        en2 = Ne
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2  ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(gaux)
                    synapses( i1, j1) = gee1 + dgee1*(-1 + 2*gaux)
                END IF
            END DO
        END DO
        !! E -> I
        st1 =  1
        en1 =  Ne
        st2 = Ne+1
        en2 = Ne+Ni
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2  ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(gaux)
                    synapses( i1, j1) = gei1 + dgei1*(-1 + 2*gaux)
                END IF
            END DO
        END DO
        !! I -> E
        st1 = Ne+1
        en1 = Ne+Ni
        st2 =  1
        en2 =  Ne
        DO i1 = st1, en1  ! Rows
            DO j1 = st2, en2 ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(gaux)
                    synapses( i1, j1) = gie1 + dgie1*(-1 + 2*gaux)
                END IF
            END DO
        END DO
        !! I -> I
        st1 =  Ne+1
        en1 =  Ne+Ni
        st2 =  Ne+1
        en2 =  Ne+Ni
        DO i1 = st1, en1
            DO j1 = st2, en2 ! Columns
                IF (A(i1,j1) .EQ. 1) THEN
                    CALL random_number(gaux)
                    synapses( i1, j1) = gii1 + dgii1*(-1 + 2*gaux)
                END IF
            END DO
        END DO
    END function synapses
    
END MODULE modu
