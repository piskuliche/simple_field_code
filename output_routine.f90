SUBROUTINE WRITE_HD5F(dot1, dot2, eoh1, eoh2, z0, nmol, nconfig)
! **********************************************************************************************************************
! This subroutine writes the HDF5 file, for all of the OHs
! 
! The hierarchy of the HDF5 file is as follows:
!   Here: a, b are indices for the first and second OHs of a particular molecule.
!
!  file: field.h5
!    - dot_a: dot product values for OH1
!    - dot_b: dot product values for OH2
!    - eoh_a: OH vector values for OH1
!    - eoh_b: OH vector values for OH2
!    - z0_a: z0 values for OH1
!    - z0_b: z0 values for OH2
! 
! Input:
!   - dot1, dot2: Dot product arrays (nmol, ntimes)
!   - eoh1, eoh2: OH vector arrays (nmol, ntimes, 3)
!   - z0: z-coordinate values (for sfg) (nmol, ntimes)
!   - nmol: number of molecules
!   - nconfig: number of configurations
! 
! Output:
!   - HDF5 file (field.h5): Contains dot products, OH vectors, and z0 values for each OH
!
! TODO:
!   1) Could maybe combine this with WRITE_HD5F_Samples by always having it write a set number of samples
!   
! **********************************************************************************************************************


    USE HDF5

    IMPLICIT NONE

    REAL, DIMENSION(:,:), INTENT(IN) :: dot1, dot2   ! Dimensions(imol, iconfig)
    REAL, DIMENSION(:,:), INTENT(IN) :: z0 ! Dimensions(imol, iconfig)
    REAL, DIMENSION(:,:,:), INTENT(IN) :: eoh1, eoh2 ! Dimensions(imol, idim, iconfig)
    INTEGER, INTENT(IN) :: nmol, nconfig

    INTEGER :: i

    ! HDF5 Variables for
    CHARACTER(LEN=8), PARAMETER :: filename = "field.h5"    ! File name
    INTEGER(HID_T) :: file_id                               ! File identifier
    INTEGER :: ERROR_FLAG! Error flag

    INTEGER(HID_T) :: dataspace_id, dataset_id
    CHARACTER(LEN=20) :: dataset_name
    INTEGER(HSIZE_T), DIMENSION(1) :: dot_dims 
    INTEGER(HSIZE_T), DIMENSION(2) :: eoh_dims 

    dot_dims = (/nconfig/)
    eoh_dims = (/3, nconfig/)

    ! Initialize FORTRAN interface

    CALL h5open_f(ERROR_FLAG)
    ! Create a new file using default properties.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ERROR_FLAG)

    ! Loop over water
    DO i=1, nmol
        ! *** HDF5 NOTE *** !
        ! * Because HDF5 is not super self-explanatory, it should be noted that the way this works
        ! * is to open the library, create the file (pre-do loop)
        ! * then create the dataspace, dataset, and write to it.
        ! * Then you close the dataspace and dataset. Then you repeat for the next dataset.
        ! * Then you close the file, and the library (post do loop)
        ! * So inside the loop, each dataset has 5 calls to the HDF5 library
        ! *** END HDF5 NOTE *** !
        
        ! Write first OH
        !   dot product value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot1(i,:), dot_dims, ERROR_FLAG)
        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   z0 value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "z0_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, z0(i,:), dot_dims, ERROR_FLAG)
        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   eoh value
        CALL h5screate_simple_f(2, eoh_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh1(i,:,:), eoh_dims, ERROR_FLAG)


        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        ! Write Second OH 
        !   dot product value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot2(i,:), dot_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   z0 value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "z0_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, z0(i,:), dot_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)
        !   eoh value
        CALL h5screate_simple_f(2, eoh_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh2(i,:,:), eoh_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)
    END DO


    ! Terminate ACCESS to the HDF5 library
    CALL h5fclose_f(file_id, ERROR_FLAG)

    ! Close Fortran interface
    CALL h5close_f(ERROR_FLAG)

END SUBROUTINE WRITE_HD5F


SUBROUTINE WRITE_HD5F_Samples(dot1, dot2, eoh1, eoh2, z0, nmol, nconfig, nsamples, samples)
! **********************************************************************************************************************
! This subroutine writes the HDF5 file, for a set number of OHs
! 
! The hierarchy of the HDF5 file is as follows:
!   Here: a, b are indices for the first and second OHs of a particular molecule.
!
!  file: field.h5
!    - dot_a: dot product values for OH1
!    - dot_b: dot product values for OH2
!    - eoh_a: OH vector values for OH1
!    - eoh_b: OH vector values for OH2
!    - z0_a: z0 values for OH1
!    - z0_b: z0 values for OH2
! 
! Input:
!   - dot1, dot2: Dot product arrays (nsamples, ntimes)
!   - eoh1, eoh2: OH vector arrays (nsamples, ntimes, 3)
!   - z0: z-coordinate values (for sfg) (nsamples, ntimes)
!   - nmol: number of molecules
!   - nconfig: number of configurations
!   - nsamples: number of samples to write
!   - samples: array of samples to write (nsamples)
! 
! Output:
!   - HDF5 file (field.h5): Contains dot products, OH vectors, and z0 values for each OH
!
! TODO:
!   1) Could maybe combine this with WRITE_HD5F by always having it write a set number of samples
!   
! **********************************************************************************************************************

    USE HDF5

    IMPLICIT NONE

    REAL, DIMENSION(:,:), INTENT(IN) :: dot1, dot2   ! Dimensions(imol, iconfig)
    REAL, DIMENSION(:,:), INTENT(IN) :: z0
    REAL, DIMENSION(:,:,:), INTENT(IN) :: eoh1, eoh2 ! Dimensions(imol, idim, iconfig)
    INTEGER, INTENT(IN) :: nmol, nconfig
    INTEGER, INTENT(IN) :: nsamples
    INTEGER, DIMENSION(:), INTENT(IN) :: samples

    INTEGER :: ival, i

    ! HDF5 Variables for
    CHARACTER(LEN=8), PARAMETER :: filename = "field.h5"    ! File name
    INTEGER(HID_T) :: file_id                               ! File identifier
    INTEGER :: ERROR_FLAG! Error flag

    INTEGER(HID_T) :: dataspace_id, dataset_id
    CHARACTER(LEN=20) :: dataset_name
    INTEGER(HSIZE_T), DIMENSION(1) :: dot_dims 
    INTEGER(HSIZE_T), DIMENSION(2) :: eoh_dims 

    dot_dims = (/nconfig/)
    eoh_dims = (/3, nconfig/)

    ! Initialize FORTRAN interface

    CALL h5open_f(ERROR_FLAG)
    ! Create a new file using default properties.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ERROR_FLAG)

    ! Loop over water
    DO i=1, nsamples
        ! *** HDF5 NOTE *** !
        ! * Because HDF5 is not super self-explanatory, it should be noted that the way this works
        ! * is to open the library, create the file (pre-do loop)
        ! * then create the dataspace, dataset, and write to it.
        ! * Then you close the dataspace and dataset. Then you repeat for the next dataset.
        ! * Then you close the file, and the library (post do loop)
        ! * So inside the loop, each dataset has 5 calls to the HDF5 library
        ! *** END HDF5 NOTE *** !

        ival = samples(i)
        ! Write first OH
        !   dot product value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot1(ival,:), dot_dims, ERROR_FLAG)
        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   dot product value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "z0_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, z0(ival,:), dot_dims, ERROR_FLAG)
        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   eoh value
        CALL h5screate_simple_f(2, eoh_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+1
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh1(ival,:,:), eoh_dims, ERROR_FLAG)


        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        ! Write Second OH 
        !   dot product value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot2(ival,:), dot_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   z0 value
        CALL h5screate_simple_f(1, dot_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "z0_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, z0(ival,:), dot_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)

        !   eoh value
        CALL h5screate_simple_f(2, eoh_dims, dataspace_id, ERROR_FLAG)
        WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+2
        CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
        CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh2(ival,:,:), eoh_dims, ERROR_FLAG)

        CALL h5dclose_f(dataset_id, ERROR_FLAG)
        CALL h5sclose_f(dataspace_id, ERROR_FLAG)
    END DO


    ! Terminate ACCESS to the HDF5 library
    CALL h5fclose_f(file_id, ERROR_FLAG)

    ! Close Fortran interface
    CALL h5close_f(ERROR_FLAG)

END SUBROUTINE WRITE_HD5F_Samples