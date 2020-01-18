!==============================================================================!
  subroutine Comm_Mod_Create_New_Types(comm)
!------------------------------------------------------------------------------!
!   Creates new data type for MPI I/O.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!
!   There is an issue with this procedure, but it's more related to MPI/IO     !
!   functions than T-Flows.  In cases a subdomain has no physical boundary     !
!   cells, variable "nb_s" turns out to be zero.  This, per se, should not     !
!   be an issue if MPI/IO functions could handle calls to:                     !
!   "Mpi_Type_Create_Indexed_Block(...)" and "Mpi_File_Write(...)" with zero   !
!   length.  But they don't.  Therefore, I avoid allocation with zero size     !
!   (max(nb_s,1)) here and creation of new types with zero size in             !
!   "Load_Maps".  It is a bit of a dirty trick :-(                             !
!------------------------------------------------------------------------------!

  ! Create new type for cells
  call Mpi_Type_Create_Indexed_Block(comm % nc_s,     & ! length of map
                                     1,               & ! block size
                                     comm % cell_map, & ! displacements 
                                     MPI_DOUBLE,      & ! old data type
                                     cell_map_type,   & ! new data type 
                                     error)             ! integer error
  call Mpi_Type_Commit(cell_map_type, error)

  ! Create new type for boundary cells
  call Mpi_Type_Create_Indexed_Block(max(comm % nb_s,1),  & ! map length
                                     1,                   & ! block size
                                     comm % bnd_cell_map, & ! displacem.
                                     MPI_DOUBLE,          & ! old type
                                     bnd_cell_map_type,   & ! new type 
                                     error)                 ! int. error
  call Mpi_Type_Commit(bnd_cell_map_type, error)

  end subroutine
