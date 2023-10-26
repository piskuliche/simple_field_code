SUBROUTINE WRITE_HD5F(dot1, dot2, eoh1, eoh2, z0, nmol, nconfig)

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