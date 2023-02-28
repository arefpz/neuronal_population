MODULE WRITER
    USE VARS
    IMPLICIT none

    CONTAINS

    subroutine write_rho_part()
        REAL :: xrand
        INTEGER :: fileID
        CHARACTER(len = 1024) :: fileaux
        CHARACTER(len=1024) :: filename1
        INTEGER(KIND = DP) :: rho_siz
        INTEGER(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: rhox
        INTEGER(KIND = DP), ALLOCATABLE, DIMENSION(:) :: rhox_ind
        INTEGER(KIND = DP) :: k2, j2, n1
        ALLOCATE(rhox(neurons, ptime),rhox_ind(neurons))
        CALL random_number(xrand)
        fileID = INT(xrand * 1000)
        rhox_ind = 0
        rhox  = 0
        rho_siz = size(rho,2)
        DO k2 = 1_DP, neurons
            n1 = 0
            DO j2 = 1_DP, rho_siz
                IF (rho(k2,j2) .EQ. 1) THEN
                    n1 = n1 + 1
                    rhox(k2,n1) = j2
                END IF
            END DO
            rhox_ind(k2) = n1
        END DO
        WRITE(fileaux,*) "(A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
        WRITE(filename1,TRIM(fileaux)) 'rho_part',INT(part),"_ens",INT(ens),".txt"     
        OPEN(UNIT=fileID, FILE=TRIM(filename1), ACTION = "WRITE", STATUS = "UNKNOWN")
        
        
        DO k2 = 1, neurons
            WRITE(fileID,*) (rhox(k2, j2) , j2 = 1, maxval(rhox_ind))
            ! WRITE(1010,*) (rhox(k2, j2) , j2 = 1, ptime)
        END DO
        CLOSE(fileID)
        PRINT*,"rho, part ", part, " is written to the disk."
    end subroutine write_rho_part


    subroutine write_A(A,neurons)
        implicit NONE
        INTEGER :: k1, J1
        INTEGER :: neurons
        INTEGER(KIND = 1), INTENT(IN) :: A(:,:)
        OPEN(UNIT = 200, FILE = "A.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
        DO k1 = 1, neurons
            WRITE(200, *) (A(k1, J1), J1 = 1, neurons)
        END DO
        CLOSE(200)
        PRINT*, "The connectivity matrix got written onto the disk."
    end subroutine write_A

    subroutine write_1DR(XSave,xname)
        IMPLICIT NONE
        REAL :: xrand
        INTEGER :: fileID
        CHARACTER(len = 1024) :: fileaux
        INTEGER :: j1
        CHARACTER(len = 1024) :: filenames
        CHARACTER(len = *), INTENT(IN) :: xname
        REAL(KIND = DP), INTENT(IN) :: XSave(:)
        CALL random_number(xrand)
        fileID = INT(xrand * 10000)    
        IF ((ens .EQ. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(0),".txt" 
        ELSEIF ((ens .EQ. 0) .AND. (part .NE. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(0),".txt" 
        ELSE IF ((ens .NE. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(ens),".txt" 
        ELSE
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(ens),".txt"
        END IF
        OPEN(UNIT = fileID, FILE = TRIM(filenames), ACTION = "WRITE", STATUS = "UNKNOWN")
        DO j1=1, size(XSave,1)
            WRITE(fileID, *) XSave(j1)
        END DO
        CLOSE(fileID)
        PRINT*, xname, " got written onto the disk."
    end subroutine write_1DR
    
    subroutine write_1DI(XSave,xname)
        IMPLICIT NONE
        REAL :: xrand
        INTEGER :: fileID
        CHARACTER(len = 1024) :: fileaux
        INTEGER :: j1
        CHARACTER(len = 1024) :: filenames
        CHARACTER(len = *), INTENT(IN) :: xname
        INTEGER(KIND = DP), INTENT(IN) :: XSave(:)
        CALL random_number(xrand)
        fileID = FLOOR(xrand * 10000)
        IF ((ens .EQ. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(0),".txt" 
        ELSEIF ((ens .EQ. 0) .AND. (part .NE. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(0),".txt" 
        ELSE IF ((ens .NE. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(ens),".txt" 
        ELSE
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(ens),".txt"
        END IF
        
        OPEN(UNIT = fileID, FILE = TRIM(filenames), ACTION = "WRITE", STATUS = "UNKNOWN")

        DO j1=1, size(XSave,1)
            WRITE(fileID, *) XSave(j1)
        END DO
        CLOSE(fileID)
        PRINT*, xname, " got written onto the disk."
    end subroutine write_1DI
    
    subroutine write_2DI(XSave,xname)
        IMPLICIT NONE
        REAL :: xrand
        INTEGER :: fileID
        CHARACTER(len = 1024) :: fileaux
        INTEGER :: j1, k1, nrow, nclo
        INTEGER(KIND = DP), INTENT(IN) :: XSave(:,:)
        CHARACTER(len = 1024) :: filenames
        CHARACTER(len = *), INTENT(IN) :: xname
        CALL random_number(xrand)
        fileID = INT(xrand * 10000)
        nrow = size(Xsave, 1)
        Nclo = size(Xsave, 2)
        IF ((ens .EQ. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(0),".txt" 
        ELSEIF ((ens .EQ. 0) .AND. (part .NE. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(0),".txt" 
        ELSE IF ((ens .NE. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(ens),".txt" 
        ELSE
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(ens),".txt"
        END IF
        OPEN(UNIT = fileID, FILE = TRIM(filenames), ACTION = "WRITE", STATUS = "UNKNOWN")

        DO j1=1, nrow
            WRITE(fileID, *) (XSave(j1,k1), k1 = 1, nclo )
        END DO
        CLOSE(fileID)
        PRINT*, xname, " got written onto the disk."
    end subroutine write_2DI
    
    subroutine write_2DR(XSave,xname)
        IMPLICIT NONE
        INTEGER :: j1, k1, nrow, nclo
        REAL :: xrand
        INTEGER :: fileID
        CHARACTER(len = 1024) :: fileaux
        REAL(KIND = 4), INTENT(IN) :: XSave(:,:)
        CHARACTER(len = 1024) :: filenames
        CHARACTER(len = *), INTENT(IN) :: xname
        CALL random_number(xrand)
        fileID = INT(xrand * 10000)
        nrow = size(Xsave, 1)
        Nclo = size(Xsave, 2)
        IF ((ens .EQ. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(0),".txt" 
        ELSEIF ((ens .EQ. 0) .AND. (part .NE. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",0,"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(0),".txt" 
        ELSE IF ((ens .NE. 0) .AND. (part .EQ. 0)) THEN
            WRITE(fileaux,*) "(A","A","I",0,&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(0),"_ens",INT(ens),".txt" 
        ELSE
            WRITE(fileaux,*) "(A","A","I",FLOOR(LOG10(REAL(part))+1),&
                         "A","I",FLOOR(LOG10(REAL(ens))+1),"A",")"
            WRITE(filenames,TRIM(fileaux)) xname,'_part',INT(part),"_ens",INT(ens),".txt"
        END IF
        OPEN(UNIT = fileID, FILE = TRIM(filenames), ACTION = "WRITE", STATUS = "UNKNOWN")
        
        DO j1=1, nrow
            WRITE(fileID, *) (XSave(j1,k1), k1 = 1, nclo )
        END DO
        CLOSE(fileID)
        PRINT*, xname, " got written onto the disk."
    end subroutine write_2DR

END MODULE WRITER
