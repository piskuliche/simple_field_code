SUBROUTINE Read_Input(nconfig, nmoltypes, molnames, nmols, natoms, charges, rmax, L, which_is_water)
    IMPLICIT NONE

    ! Output Variables
    REAL, INTENT(OUT) :: rmax  ! Cutoff radius (in Angstroms)
    REAL, DIMENSION(3), INTENT(OUT) :: L ! Box Length (in Angstroms)

    INTEGER :: i
    INTEGER :: nmoltypes, nconfig, which_is_water, nwater
    INTEGER, DIMENSION(10) :: nmols, natoms
    REAL, DIMENSION(10, 2000) :: charges
    CHARACTER(LEN=10), DIMENSION(10) :: molnames
    
    ! Read the field input file
    OPEN(10, file='field_input.in', status='old')
        READ(10,*)
        READ(10,*) nmoltypes, nconfig, which_is_water
        READ(10,*) 
        READ(10,*) rmax, L(1), L(2), L(3)
        READ(10,*)
        IF (nmoltypes > 10) STOP 'Too many molecule types'
        READ(10,*) (molnames(i), i=1, nmoltypes) ! molecule names
        READ(10,*)
        READ(10,*) (nmols(i), i=1, nmoltypes) ! number of mols
    CLOSE(10)

    WRITE(*,*) "nm1", nmols(1)

    DO i=1, nmoltypes
        CALL Read_Molecule(i, molnames(i), charges(i,:), natoms(i))
        IF ( i == which_is_water) THEN
            nwater = nmols(i)
        END IF
    ENDDO

END SUBROUTINE Read_Input

SUBROUTINE Read_Molecule(imol, molname, q, natoms)
    IMPLICIT NONE

    CHARACTER(LEN=10), INTENT(IN) :: molname
    INTEGER, INTENT(IN) :: imol
    INTEGER, INTENT(OUT) :: natoms
    REAL, DIMENSION(2000) :: q

    INTEGER :: i, ifile

    ifile = imol + 10

    OPEN(ifile, file=trim(molname)//".in", status='old')
        READ(ifile,*) natoms
        DO i=1, natoms
            READ(ifile,*) q(i)
        END DO
    CLOSE(ifile)

END SUBROUTINE Read_Molecule

SUBROUTINE Open_Field_Files(nwater)
! This subroutine opens the nwater field files
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nwater

    INTEGER ifile, i
    CHARACTER(LEN=4), DIMENSION(2*nwater) :: ext

    ! Field File extensions
    DO i =1, 2*nwater
        write(ext(i), '(I0.4)') i+100
    ENDDO

    ! Open Field Files
    DO ifile = 101, 2*nwater+100
        OPEN(ifile,file='field.'//ext(ifile-100))
    ENDDO

END SUBROUTINE Open_Field_Files

SUBROUTINE Close_Field_Files(nwater)
! This subroutine closes the nwater field files
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nwater

    INTEGER :: ifile

    DO ifile=101, 2*nwater+100
        CLOSE(ifile)
    ENDDO

END SUBROUTINE Close_Field_Files

SUBROUTINE Read_Trajectory(nconfig, nmoltypes, nmols, natoms, which_is_wat, L, rO, r1, r2, rmol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nconfig, nmoltypes, nmols(:), natoms(:), which_is_wat
    REAL, DIMENSION(3), INTENT(IN) :: L
    REAL, INTENT(INOUT) :: rO(:,:,:), r1(:,:,:), r2(:,:,:), rmol(:,:,:,:,:)
    INTEGER :: i, k, z, type, jatom
    CHARACTER(LEN=10) :: ctmp

    OPEN(11, file='../traj.xyz', status='old')

    DO z=1, nconfig
        DO i=1,2
            read(11,*) 
        END DO

        DO type=1, nmoltypes
            WRITE(*,*) "type", type
            WRITE(*,*) "nmols", nmols(type)
            IF (type == which_is_wat) THEN
                DO i=1, nmols(type)
                    ! Reading rule for water
                    read(11,*) ctmp, (rO(i,k,z), k=1,3)
                    read(11,*) ctmp, (r1(i,k,z), k=1,3)
                    read(11,*) ctmp, (r2(i,k,z), k=1,3)
                    ! Make the molecule whole
                    DO k=1, 3
                        r1(i,k,z) = r1(i,k,z) - L(k)*anint((r1(i,k,z)-rO(i,k,z))/L(k))
                        r2(i,k,z) = r2(i,k,z) - L(k)*anint((r2(i,k,z)-rO(i,k,z))/L(k))
                    ENDDO ! k
                ENDDO ! i
            ELSE ! (type == which_is_wat)
                ! Reading rule for not water
                DO i=1, nmols(type)
                    DO jatom=1, natoms(i)
                        read(11,*) ctmp, (rmol(type,i,jatom,k,z), k=1,3)
                        ! Make the molecule whole
                        IF (jatom > 1) THEN
                            DO k=1, 3
                                rmol(type,i,jatom,k,z) = rmol(type,i,jatom,k,z) &
                                & - L(k)*anint((rmol(type,i,jatom,k,z)-rmol(type,i,1,k,z))/L(k))
                            ENDDO ! k
                        END IF ! (jatom > 1)
                    END DO !jatom
                ENDDO ! i
            END IF ! type == which_is_wat
        END DO ! type
    END DO ! z -> end read_traj
    CLOSE(11)

END SUBROUTINE Read_Trajectory