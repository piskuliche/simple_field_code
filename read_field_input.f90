SUBROUTINE Read_Input(nconfig, nmoltypes, molnames, nmols, natoms, charges, rmax, L, which_is_water, nsamples)
    IMPLICIT NONE

    ! Output Variables
    REAL, INTENT(OUT) :: rmax  ! Cutoff radius (in Angstroms)
    REAL, DIMENSION(3), INTENT(OUT) :: L ! Box Length (in Angstroms)

    INTEGER :: i
    INTEGER :: nmoltypes, nconfig, which_is_water, nwater
    INTEGER, DIMENSION(10), INTENT(INOUT) :: nmols, natoms
    INTEGER, INTENT(OUT) :: nsamples
    REAL, DIMENSION(10, 2000) :: charges
    CHARACTER(LEN=10), DIMENSION(10) :: molnames

    nmols = 0; natoms = 0
    
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
        READ(10,*) 
        READ(10,*) nsamples
    CLOSE(10)

    IF (nsamples > nmols(which_is_water)) THEN
        WRITE(*,*) 'Error: nsamples > nwater'
        STOP 'Too many samples'
    END IF

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
    CHARACTER(LEN=5), DIMENSION(2*nwater) :: ext

    ! Field File extensions
    DO i =1, 2*nwater
        write(ext(i), '(I0.5)') i+100
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
    INTEGER, INTENT(IN) :: nconfig, nmoltypes, which_is_wat
    INTEGER, DIMENSION(10), INTENT(IN) :: nmols, natoms
    REAL, DIMENSION(3), INTENT(IN) :: L
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: rO, r1, r2
    REAL, DIMENSION(:,:,:,:,:), INTENT(INOUT) :: rmol
    INTEGER :: i, k, z, type, jatom
    CHARACTER(LEN=10) :: ctmp

    OPEN(11, file='../traj.xyz', status='old')

    DO z=1, nconfig
        DO i=1,2
            read(11,*) 
        END DO

        DO type=1, nmoltypes
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
                    DO jatom=1, natoms(type)
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

SUBROUTINE Read_Samples(nsamples, samples)
! ************************************************************************************
! This subroutine reads the samples.in file
! It contains the indices of the water molecules to be used for the calculation
! This is a cost saving tool for dealing with the expense of very big systems.
! In:
!  -nsamples: number of samples
! Out:
!  -samples: array of indices of the water molecules to be used for the calculation
! ************************************************************************************
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nsamples
    INTEGER, DIMENSION(:), INTENT(OUT) :: samples

    INTEGER :: i

    OPEN(15, file="samples.in", status='old')

    DO i=1, nsamples
        READ(15,*) samples(i)
    END DO

    CLOSE(15)

END SUBROUTINE Read_Samples

SUBROUTINE Read_Trajectory_Frames(unit, ntoread, nmoltypes, nmols, natoms, which_is_wat, L, rO, r1, r2, rmol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit, ntoread, nmoltypes, which_is_wat
    INTEGER, DIMENSION(10), INTENT(IN) :: nmols, natoms
    REAL, DIMENSION(3), INTENT(IN) :: L
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: rO, r1, r2
    REAL, DIMENSION(:,:,:,:,:), INTENT(INOUT) :: rmol
    INTEGER :: i, k, z, type, jatom
    CHARACTER(LEN=10) :: ctmp


    DO z=1, ntoread
        CALL Read_XYZ_Frame(unit, nmoltypes, nmols, natoms, which_is_wat, L, rO(:,:,z), r1(:,:,z), r2(:,:,z), rmol(:,:,:,:,z))
    END DO ! z -> end read_traj
    CLOSE(11)

END SUBROUTINE Read_Trajectory

SUBROUTINE Read_XYZ_Frame(unit, nmoltypes, nmols, natoms, which_is_wat, L, rO, r1, r2, rmol)
! ************************************************************************************
! This subroutine reads a frame from the trajectory file in the XYZ format
! In:
!  -nmoltypes: number of molecule types
!  -nmols: array of number of molecules of each type
!  -natoms: array of number of atoms of each type
!  -which_is_wat: index of the water molecule type
!  -L: box length
! Out:
!  -rO: array of oxygen positions
!  -r1: array of hydrogen 1 positions
!  -r2: array of hydrogen 2 positions
!  -rmol: array of all the other atoms positions
! ************************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(in) :: nmoltypes
    INTEGER, DIMENSION(:), INTENT(IN) :: nmols, natoms
    INTEGER :: which_is_wat
    REAL, DIMENSION(3) :: L
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: rO, r1, r2
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: rmol

    INTEGER :: i, type
    CHARACTER(LEN=10) :: ctmp

    rO = 0.0; r1 = 0.0; r2 = 0.0
    rmol = 0.0

    DO i=1,2
        READ(unit,*) 
    END DO

    DO type=1, nmoltypes
        IF (type == which_is_wat) THEN
            DO i=1, nmols(type)
                ! Reading rule for water
                read(unit,*) ctmp, (rO(i,k), k=1,3)
                read(unit,*) ctmp, (r1(i,k), k=1,3)
                read(unit,*) ctmp, (r2(i,k), k=1,3)
                ! Make the molecule whole
                DO k=1, 3
                    r1(i,k) = r1(i,k,z) - L(k)*anint((r1(i,k)-rO(i,k))/L(k))
                    r2(i,k) = r2(i,k,z) - L(k)*anint((r2(i,k)-rO(i,k))/L(k))
                ENDDO ! k
            ENDDO ! i
        ELSE ! (type == which_is_wat)
            ! Reading rule for not water
            DO i=1, nmols(type)
                DO jatom=1, natoms(type)
                    read(unit,*) ctmp, (rmol(type,i,jatom,k), k=1,3)
                    ! Make the molecule whole
                    IF (jatom > 1) THEN
                        DO k=1, 3
                            rmol(type,i,jatom,k,z) = rmol(type,i,jatom,k) &
                            & - L(k)*anint((rmol(type,i,jatom,k)-rmol(type,i,1,k))/L(k))
                        ENDDO ! k
                    END IF ! (jatom > 1)
                END DO !jatom
            ENDDO ! i
        END IF ! type == which_is_wat
    END DO ! type
END SUBROUTINE