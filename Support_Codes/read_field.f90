SUBROUTINE read_field(ioh)
    USE HDF5

    IMPLICIT NONE
    INTEGER :: ioh, j, k
    REAL, DIMENSION(2000) :: etmp
    DOUBLE PRECISION ::  muprime, x01tmp, x12tmp
    !DOUBLE PRECISION, DIMENSION(ntimes) :: w01, w12
    !DOUBLE PRECISION, DIMENSION(ntimes) :: mu01, mu12
    !DOUBLE PRECISION, DIMENSION(ntimes,3) :: eoh
    
    ! HDF5 Variables
    CHARACTER(len=8), PARAMETER :: filename = 'field.h5'
    INTEGER(HID_T) :: file_id, dataset_id
    INTEGER :: ERROR_FLAG

    INTEGER(HSIZE_T), DIMENSION(1) :: dot_dims 
    INTEGER(HSIZE_T), DIMENSION(2) :: eoh_dims 


    CHARACTER(len=20) :: dataset_name

    ! Open the hdf5 library
    CALL h5open_f(ERROR_FLAG)

    ! Open the hdf5 file
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ERROR_FLAG)

    ! Set dataset name
    WRITE(dataset_name, '(A, I0)') 'dot_', ioh
    
    ! Open the existing dataset
    CALL h5dopen_f(file_id, dataset_name, dataset_id, ERROR_FLAG)

    ! Read the dataset
    CALL h5dread_f(dataset_id, H5T_NATIVE_REAL, etmp, dot_dims, ERROR_FLAG)

    ! Close the dataset
    CALL h5dclose_f(dataset_id, ERROR_FLAG)

    ! Close the file
    CALL h5fclose_f(file_id, ERROR_FLAG)
    
    ! Close the library
    CALL h5close_f(ERROR_FLAG)

    WRITE(*,*) "Read value", etmp(1), etmp(2)

END SUBROUTINE read_field