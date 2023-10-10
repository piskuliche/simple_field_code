SUBROUTINE WRITE_HD5F(dot1, dot2, eoh1, eoh2, nmol, nconfig)

USE HDF5

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: dot1, dot2   ! Dimensions(imol, iconfig)
REAL, DIMENSION(:,:,:), INTENT(IN) :: eoh1, eoh2 ! Dimensions(imol, idim, iconfig)
INTEGER, INTENT(IN) :: nmol, nconfig

INTEGER :: i

! HDF5 Variables for
CHARACTER(LEN=8), PARAMETER :: filename = "field.h5"    ! File name
INTEGER(HID_T) :: file_id                               ! File identifier
INTEGER :: ERROR_FLAG! Error flag

INTEGER :: dataspace_id, dataset_id
CHARACTER(LEN=20) :: dataset_name

! Initialize FORTRAN interface

CALL h5open_f(ERROR_FLAG)
! Create a new file using default properties.
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ERROR_FLAG)

! Loop over water
DO i=1, nmol
    ! Write first OH
    !   dot product value
    CALL h5screate_simple(1, [nconfig], dataspace_id, ERROR_FLAG)
    WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+1
    CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot1(i,:), ERROR_FLAG)

    CALL h5dclose_f(dataset_id)
    CALL h5sclose_f(dataspace_id)

    !   eoh value
    CALL h5screate_simple(1, [3, nconfig], dataspace_id, ERROR_FLAG)
    WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+1
    CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh1(i,:,:), ERROR_FLAG)

    CALL h5dclose_f(dataset_id)
    CALL h5sclose_f(dataspace_id)

    ! Write Second OH 
    !   dot product value
    CALL h5screate_simple(1, [nconfig], dataspace_id, ERROR_FLAG)
    WRITE(dataset_name, '(A, I0)') "dot_", (i-1)*2+2
    CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dot2(i,:), ERROR_FLAG)

    CALL h5dclose_f(dataset_id)
    CALL h5sclose_f(dataspace_id)
    !   eoh value
    CALL h5screate_simple(1, [3, nconfig], dataspace_id, ERROR_FLAG)
    WRITE(dataset_name, '(A, I0)') "eoh_", (i-1)*2+1
    CALL h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dataspace_id, dataset_id, ERROR_FLAG)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, eoh2(i,:,:), ERROR_FLAG)

    CALL h5dclose_f(dataset_id)
    CALL h5sclose_f(dataspace_id)
END DO


! Terminate ACCESS to the HDF5 library
CALL h5fclose_f(file_id, ERROR_FLAG)

! Close Fortran interface
CALL h5close_f(ERROR_FLAG)

END SUBROUTINE WRITE_HD5F