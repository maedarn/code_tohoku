PROGRAM TEST
  USE HDF5 ! This module contains all necessary HDF5 modules
  IMPLICIT NONE

  ! Names (file and HDF5 objects)
  CHARACTER(LEN=15), PARAMETER :: filename = "mytestfile.hdf5" ! File name
  CHARACTER(LEN=9), PARAMETER :: dsetname1 = "mydataset" ! Dataset name
  CHARACTER(LEN=8), PARAMETER :: groupname = "subgroup" ! Sub-Group 1 name
  CHARACTER(LEN=9), PARAMETER :: groupname3 = "subgroup2" ! Sub-Group 3 name
  ! Dataset 2 name
  CHARACTER(LEN=24), PARAMETER :: dsetname2 = "subgroup/another_dataset"
  ! Dataset 3 name
  CHARACTER(LEN=23), PARAMETER :: dsetname3 = "subgroup2/dataset_three"

  ! Identifiers
  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: group_id      ! Group identifier
  INTEGER(HID_T) :: group3_id     ! Group 3 identifier
  INTEGER(HID_T) :: dset1_id      ! Dataset 1 identifier
  INTEGER(HID_T) :: dset2_id      ! Dataset 2 identifier
  INTEGER(HID_T) :: dset3_id      ! Dataset 3 identifier
  INTEGER(HID_T) :: dspace1_id    ! Dataspace 1 identifier
  INTEGER(HID_T) :: dspace2_id    ! Dataspace 2 identifier
  INTEGER(HID_T) :: dspace3_id    ! Dataspace 3 identifier

  ! Integer array
  INTEGER :: rank                 ! Dataset rank
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1 = (/100/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
  INTEGER, DIMENSION(100) :: dset_data1   ! Data buffers

  ! FP array
  INTEGER(HSIZE_T), DIMENSION(1) :: dims2 = (/50/)
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2
  REAL, DIMENSION(50) :: dset_data2

  ! Integer array for dataset_three
  INTEGER(HSIZE_T), DIMENSION(1) :: dims3 = (/10/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims3     ! Dataset rank
  INTEGER, DIMENSION(10) :: dset_data3

  ! Misc variables (e.g. loop counters)
  INTEGER :: error ! Error flag
  INTEGER :: i,j
  ! =====================================================================

  ! Initialize the dset_data array 
  data_dims1(1) = 100
  rank = 1
  DO i = 1, 100
     dset_data1(i) = i
  END DO

  ! Initialize Fortran interface
  CALL h5open_f(error)   
  ! Create a new file
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

  ! Create dataspace 1 (the dataset is next) "dspace_id" is returned
  CALL h5screate_simple_f(rank, dims1, dspace1_id, error)
  ! Create dataset 1 with default properties "dset_id" is returned
  CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_INTEGER, dspace1_id, &
       dset1_id, error)
  ! Write dataset 1
  CALL h5dwrite_f(dset1_id, H5T_NATIVE_INTEGER, dset_data1, data_dims1, &
       error)
  ! Close access to dataset 1
  CALL h5dclose_f(dset1_id, error)
  ! Close access to data space 1
  CALL h5sclose_f(dspace1_id, error)

  ! Create a group in the HDF5 file
  CALL h5gcreate_f(file_id, groupname, group_id, error)
  ! Close the group
  CALL h5gclose_f(group_id, error)

  ! Create dataspace 2 (the dataset is next)
  data_dims2(1) = 50
  DO i = 1, 50
     dset_data2(i) = 1.0
  END DO
  ! Create dataspace 2
  CALL h5screate_simple_f(rank, dims2, dspace2_id, error)
  ! Create dataset 2 with default properties
  CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_REAL, dspace2_id, &
       dset2_id, error)
  ! Write dataset 2
  CALL h5dwrite_f(dset2_id, H5T_NATIVE_REAL, dset_data2, data_dims2, &
       error)
  ! Close access to dataset 2
  CALL h5dclose_f(dset2_id, error)
  ! Close access to data space 2
  CALL h5sclose_f(dspace2_id, error)

  ! Create a group in the HDF5 file
  CALL h5gcreate_f(file_id, groupname3, group3_id, error)
  ! Close the group
  CALL h5gclose_f(group3_id, error)

  ! Create dataspace 3
  data_dims3(1) = 10
  DO i = 1, 10
     dset_data3(i) = i + 3
  END DO
  ! Create dataspace 3
  CALL h5screate_simple_f(rank, dims3, dspace3_id, error)
  ! Create dataset 3 with default properties
  CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_INTEGER,  &
       dspace3_id, dset3_id, error)
  ! Write dataset 3
  CALL h5dwrite_f(dset3_id, H5T_NATIVE_INTEGER, dset_data3, data_dims3, &
       error)
  ! Close access to dataset 3
  CALL h5dclose_f(dset3_id, error)
  ! Close access to data space 3
  CALL h5sclose_f(dspace3_id, error)

  ! Close the file
  CALL h5fclose_f(file_id, error)
  ! Close FORTRAN interface
  CALL h5close_f(error)
END PROGRAM TEST
