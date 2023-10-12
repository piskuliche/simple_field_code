PROGRAM Field
    use field_module
    IMPLICIT None

    ! Variables
    INTEGER :: i, imol, z            ! Loop index

    ! Input Variables
    INTEGER :: nconfig       ! Number of configurations
    INTEGER :: nmoltypes     ! Number of molecule types
    INTEGER, DIMENSION(10) :: nmols, natoms ! number of mols, and number of atoms/mol
    INTEGER :: which_is_wat  ! Which molecule is water    
    REAL, DIMENSION(10, 2000) :: charges ! Charges
    CHARACTER(LEN=10), DIMENSION(10) :: molnames ! Molecule names
    REAL :: rmax            ! Cutoff radius
    REAL, DIMENSION(3) :: L ! Box length
    INTEGER :: nsamples
    INTEGER, ALLOCATABLE, DIMENSION(:) :: samples ! Samples for field calculation
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rO, r1, r2 ! Water coordinate arrays
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: rmol ! Molecular coordinate arrays

    ! Reading Variables
    INTEGER :: max_mol, max_natom

    ! Timing Variables
    REAL :: tmp1, tmp2

    ! Field Variables
    REAL, ALLOCATABLE, DIMENSION(:,:) :: dot1, dot2
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: eOH1, eOH2

    ! Openmp
    INTEGER :: max_omp_threads, threadno

! I. Setup *******************************************************************

    CALL cpu_time(tmp1)

    ! I.A: Read the Input File ****************************************************
        ! Read the input file and grab the data that describes both the system, and the
        ! field that we are going to calculate. 
        ! Note: for each molecule, need a foo.in file, where "foo" is the molecule name.
        ! This opens all the field files as well - close them at the end.
    CALL Read_Input(nconfig, nmoltypes, molnames, nmols, natoms, charges, rmax, L, which_is_wat, nsamples)
    write(*,*) "There are ", nmoltypes, " molecule types"
    WRITE(*,*) "There are ", nconfig, " configurations"

    ! Check if using sampling
    IF (nsamples > 0) THEN 
        ALLOCATE(samples(nsamples))
        CALL Read_Samples(nsamples, samples)
    END IF


    max_mol = 0; max_natom = 0

    ! I.B: Allocate Arrays ********************************************************
        ! Allocate the arrays that will hold the coordinates of the molecules
    DO i = 1, nmoltypes
        WRITE(*,*) "There are ", nmols(i), " molecules of type ", molnames(i)
        IF (i .eq. which_is_wat) THEN
            ALLOCATE(rO(nmols(i), 3, nconfig))
            ALLOCATE(r1(nmols(i), 3, nconfig))
            ALLOCATE(r2(nmols(i), 3, nconfig))
        ELSE
            max_mol = MAX(max_mol, nmols(i))
            max_natom = MAX(max_natom, natoms(i))
        END IF
    END DO
    ALLOCATE(rmol(nmoltypes, max_mol, max_natom, 3, nconfig))
    ! Zero trajectory arrays
    rO=0d0; r1=0d0; r2=0d0; rmol=0d0

    ! Allocate and zero output arrays
    ALLOCATE(dot1(nmols(which_is_wat),nconfig))
    ALLOCATE(dot2(nmols(which_is_wat),nconfig))
    ALLOCATE(eOH1(nmols(which_is_wat),3,nconfig))
    ALLOCATE(eOH2(nmols(which_is_wat),3,nconfig))
    dot1=0d0; dot2=0d0; eOH1=0d0; eOH2=0d0

    WRITE(*,*) "wat", nmols(which_is_wat)

    ! I.C: Setup Field Files *****************************************************
        ! Open the field files and write the header
    !CALL Open_Field_Files(nmols(which_is_wat))

    CALL cpu_time(tmp2)
    WRITE(6,'(A,F10.2,A)') ' Setup time = ',tmp2-tmp1,' s'

! II. Trajectory Reading ******************************************************
    CALL cpu_time(tmp1)

    ! Opens the trajectory file and reads the coordinates of the molecules
    ! It also makes the molecules whole along the way.
    Call Read_Trajectory(nconfig, nmoltypes, nmols, natoms, which_is_wat, L, rO, r1, r2, rmol) 

    CALL cpu_time(tmp2)
    WRITE(6,'(A,F10.2,A)') ' Read time = ',tmp2-tmp1,' s'
    CALL flush(6)

! III. Field Calculation *******************************************************

    WRITE(6,*) ' Beginning field calc'
    CALL cpu_time(tmp1)
    CALL flush(6)

    ! Check if using sampling or not. Eventually just build samples to be its own array.
    IF (nsamples > 0) THEN
        Call Get_Field_Samples(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, samples, &
            & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    ELSE
        Call Get_Field(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, &
            & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    END IF

    write(*,*) dot1(1,1),"t"
    
    CALL cpu_time(tmp2)
    WRITE(6,'(A,F10.2,A)') ' Efield time = ',tmp2-tmp1,' s'


! IV. Write Data ************************************************************

    WRITE(6,*) ' nwat = ', nmols(which_is_wat)
    write(*,*) dot1(1,1),"t"
    IF (nsamples > 0) THEN
        CALL WRITE_HD5F_Samples(dot1, dot2, eoh1, eoh2, nmols(which_is_wat), nconfig, nsamples, samples)
    ELSE
        CALL WRITE_HD5F(dot1, dot2, eoh1, eoh2, nmols(which_is_wat), nconfig)
    ENDIF 


! V. Deallocate ***************************************************************

    DEALLOCATE(rO, r1, r2, rmol)
    DEALLOCATE(dot1, dot2, eOH1, eOH2)

    !CALL Close_Field_Files(nmols(which_is_wat))
    
END PROGRAM FIELD