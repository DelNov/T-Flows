!==============================================================================!
  subroutine Grid_Mod_Form_Maps(grid)
!------------------------------------------------------------------------------!
!   Forms maps for parallel backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt
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

  ! Initialize number of cells in subdomain
  grid % comm % nc_s = grid % n_cells - grid % comm % n_buff_cells

  ! Initialize number of boundary cells in subdomain
  do c = -grid % n_bnd_cells, -1
    if(grid % comm % cell_proc(c) .eq. this_proc) then
      grid % comm % nb_s = grid % comm % nb_s + 1
    end if
  end do

  ! First and last cell to send
  grid % comm % nb_f = grid % n_bnd_cells
  grid % comm % nb_l = grid % n_bnd_cells - grid % comm % nb_s + 1

  ! Initialize total number of cells
  grid % comm % nc_t = grid % comm % nc_s
  grid % comm % nb_t = grid % comm % nb_s

  !--------------------------------!
  !                                !
  !   For run with one processor   !
  !                                !
  !--------------------------------!
  if(n_proc < 2) then

    !-------------------------------------!
    !   Global cell numbers for T-Flows   !
    !-------------------------------------!
    do c = -grid % comm % nb_t, grid % comm % nc_t
      grid % comm % cell_glo(c) = c
    end do

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!
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

  !-----------------------!
  !                       !
  !   For parallel runs   !
  !                       !
  !-----------------------!
  else

    call Comm_Mod_Global_Sum_Int(grid % comm % nc_t)
    call Comm_Mod_Global_Sum_Int(grid % comm % nb_t)

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!
    allocate(grid % comm % cell_map    (grid % comm % nc_s))
    allocate(grid % comm % bnd_cell_map(max(grid % comm % nb_s,1)))
    grid % comm % cell_map(:)     = 0
    grid % comm % bnd_cell_map(:) = 0

    !---------------------!
    !   Inside cell map   !
    !---------------------!
    do c = 1, grid % comm % nc_s
      ! Take cell mapping to be the same as global cell numbers but start from 0
      grid % comm % cell_map(c) = grid % comm % cell_glo(c) - 1
    end do

    !-----------------------!
    !   Boundary cell map   !
    !-----------------------!
    cnt = 0
    do c = -grid % comm % nb_f, -grid % comm % nb_l
      cnt = cnt + 1
      grid % comm % bnd_cell_map(cnt) = grid % comm % cell_glo(c)  &
                                      + grid % comm % nb_t
    end do
  end if

  end subroutine
