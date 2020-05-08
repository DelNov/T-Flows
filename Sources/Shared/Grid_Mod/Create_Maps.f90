!==============================================================================!
  subroutine Grid_Mod_Create_Maps(grid)
!------------------------------------------------------------------------------!
!   Reads: name.map file                                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!
!   There is an issue with this procedure, but it's more related to MPI/IO     !
!   functions than T-Flows.  In cases a subdomain has no physical boundary     !
!   cells, variable "nb_s" turns out to be zero.  This, per se, should not     !
!   be an issue if MPI/IO functions could handle calls to:                     !
!   "Mpi_Type_Create_Indexed_Block(...)" and "Mpi_File_Write(...)" with zero   !
!   length.  But they don't.  Therefore, I avoid allocation with zero size     !
!   (max(nb_s,1)) here and creation of new types with zero size in             !
!   "Comm_Mod_Create_New_Types".  It is a bit of a dirty trick :-(             !
!------------------------------------------------------------------------------!

  !------------------------------------------------------------------------!
  !                                                                        !
  !   For run with one processor, no needd to read the map, just form it   !
  !                                                                        !
  !------------------------------------------------------------------------!
  if(n_proc < 2) then

    grid % comm % nc_s = grid % n_cells
    grid % comm % nb_s = grid % n_bnd_cells
    grid % comm % nc_t = grid % comm % nc_s
    grid % comm % nb_t = grid % comm % nb_s

    ! Fill up global cell numbers
    do c = -grid % comm % nb_t, grid % comm % nc_t
      grid % comm % cell_glo(c) = c
    end do

    !-----------------------------------------!
    !   Global cell numbers for MPI mapping   !
    !-----------------------------------------!
    allocate(grid % comm % cell_map    (grid % comm % nc_s))
    allocate(grid % comm % bnd_cell_map(grid % comm % nb_s))

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, grid % comm % nc_t
      grid % comm % cell_map(c) = c - 1
    end do

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, grid % comm % nb_t
      grid % comm % bnd_cell_map(c) = c - 1
    end do

  !-------------------------------------------------!
  !                                                 !
  !   For parallel runs, you need to read the map   !
  !                                                 !
  !-------------------------------------------------!
  else

    ! Map sizes
    grid % comm % nc_s = grid % n_cells
    grid % comm % nb_s = grid % n_bnd_cells

    grid % comm % nc_t  = grid % comm % nc_s
    grid % comm % nb_t  = grid % comm % nb_s
    call Comm_Mod_Global_Sum_Int(grid % comm % nc_t)
    call Comm_Mod_Global_Sum_Int(grid % comm % nb_t)

    !-----------------------------------------!
    !   Global cell numbers for MPI mapping   !
    !-----------------------------------------!
    allocate(grid % comm % cell_map    (grid % comm % nc_s))   
    allocate(grid % comm % bnd_cell_map(max(grid % comm % nb_s,1)))
    grid % comm % cell_map(:)     = 0
    grid % comm % bnd_cell_map(:) = 0

    !------------------------------------------!
    !   Take cell mapping to be the same as    !
    !   global cell numbers but start from 0   !
    !------------------------------------------!
    do c = 1, grid % comm % nc_s
      grid % comm % cell_map(c) = grid % comm % cell_glo(c) - 1
    end do

    !----------------------------------------------------------------!
    !   Correct boundary cell mapping.                               !
    !   - First it is in positive range, so insted of -nb_s to -1,   !
    !     it goes from 1 to nb_s.  (Therefore the "c+nb_s+1")        !
    !   - Second, mapping must be positive and start from zero.      !
    !    (The "+ nb_t")                                              !
    !----------------------------------------------------------------!
    do c = -grid % comm % nb_s, -1
      grid % comm % bnd_cell_map(c+grid % comm % nb_s+1) =  &
        grid % comm % cell_glo(c) + grid % comm % nb_t
    end do

    !---------------------------------------------!
    !   Refresh buffers for global cell numbers   !
    !---------------------------------------------!
    call Grid_Mod_Exchange_Int(grid, grid % comm % cell_glo)

  end if

  end subroutine
