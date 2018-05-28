!==============================================================================!
  subroutine Comm_Mod_Create_New_Types()
!------------------------------------------------------------------------------!
!   Creates new data type for MPI I/O.                                         !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Create new type for cells
  call Mpi_Type_Create_Indexed_Block(nc_s,           &  ! length of map
                                     1,              &  ! size of the block
                                     cell_map,       &  ! array of displacements 
                                     MPI_DOUBLE,     &  ! old data type
                                     cell_map_type,  &  ! new data type 
                                     error)             ! integer error
  call Mpi_Type_Commit(cell_map_type, error)

  ! Create new type for boundary cells
  call Mpi_Type_Create_Indexed_Block(nb_s,               &  ! length of map
                                     1,                  &  ! size of the block
                                     bnd_cell_map,       &  ! array of displac.
                                     MPI_DOUBLE,         &  ! old data type
                                     bnd_cell_map_type,  &  ! new data type 
                                     error)                 ! integer error
  call Mpi_Type_Commit(bnd_cell_map_type, error)

  ! Create new type for faces
  call Mpi_Type_Create_Indexed_Block(nf_s,           &  ! length of map
                                     1,              &  ! size of the block
                                     face_map,       &  ! array of displacements 
                                     MPI_DOUBLE,     &  ! old data type
                                     face_map_type,  &  ! new data type 
                                     error)             ! integer error
  call Mpi_Type_Commit(face_map_type, error)

  ! Create new type for faces
  call Mpi_Type_Create_Indexed_Block(nbf_s,              &  ! length of map
                                     1,                  &  ! size of the block
                                     buf_face_map,       &  ! array of displac.
                                     MPI_DOUBLE,         &  ! old data type
                                     buf_face_map_type,  &  ! new data type 
                                     error)                 ! integer error
  call Mpi_Type_Commit(buf_face_map_type, error)

  end subroutine
