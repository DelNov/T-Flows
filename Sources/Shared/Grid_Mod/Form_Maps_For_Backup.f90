!==============================================================================!
  subroutine Form_Maps_For_Backup(Grid)
!------------------------------------------------------------------------------!
!   Forms maps for parallel backup                                             !
!                                                                              !
!   There is an issue with this procedure, but it's more related to MPI/IO     !
!   functions than T-Flows.  In cases a subdomain has no physical boundary     !
!   cells, variable "nb_sub" turns out to be zero.  This, per se, should not   !
!   be an issue if MPI/IO functions could handle calls to:                     !
!   "Mpi_Type_Create_Indexed_Block(...)" and "Mpi_File_Write(...)" with zero   !
!   length.  But they don't.  Therefore, I avoid allocation with zero size     !
!   (max(nb_sub,1)) here and creation of new types with zero size in           !
!   "Comm_Mod_Create_New_Types".  It is a bit of a dirty trick :-(             !
!                                                                              !
!   Another issue, also related to MPI, is that maps for saving must be in     !
!   increasing order.  If you are doing only MPI that is fine, but if you      !
!   are renumbering cells for threading with OpenMP, things get cumbersome.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt, reg
!==============================================================================!

  ! Initialize number of cells in subdomain
  Grid % Comm % nc_sub = Grid % n_cells - Grid % Comm % n_buff_cells

  ! Initialize number of boundary cells in subdomain
  Grid % Comm % nb_sub = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      Grid % Comm % nb_sub = Grid % Comm % nb_sub + 1
    end if
  end do

  ! First and last cell to send are found through regions
  Grid % Comm % nb_l = Grid % n_bnd_cells
  Grid % Comm % nb_f = 1
  do reg = Boundary_Regions()
    if(Grid % region % l_cell(reg) >= Grid % region % f_cell(reg)) then
      Grid % Comm % nb_l = min(Grid % Comm % nb_l, -Grid % region % l_cell(reg))
      Grid % Comm % nb_f = max(Grid % Comm % nb_f, -Grid % region % f_cell(reg))
    end if
  end do

  ! Initialize total number of cells
  Grid % Comm % nc_tot = Grid % Comm % nc_sub
  Grid % Comm % nb_tot = Grid % Comm % nb_sub
  call Comm_Mod_Global_Sum_Int(Grid % Comm % nc_tot)
  call Comm_Mod_Global_Sum_Int(Grid % Comm % nb_tot)

  ! Allocate memory for mapping matrices
  allocate(Grid % Comm % cell_map        (Grid % Comm % nc_sub))
  allocate(Grid % Comm % bnd_cell_map(max(Grid % Comm % nb_sub,1)))
  Grid % Comm % cell_map(:)         = 0
  Grid % Comm % bnd_cell_map(:)     = 0

  !--------------------------------!
  !                                !
  !   For run with one processor   !
  !                                !
  !--------------------------------!
  if(n_proc < 2) then

    !-------------------------------------!
    !   Global cell numbers for T-Flows   !
    !-------------------------------------!
    do c = -Grid % Comm % nb_tot, Grid % Comm % nc_tot
      Grid % Comm % cell_glo(c) = c
    end do

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, Grid % Comm % nc_tot
      Grid % Comm % cell_map(c) = int(c-1, SP)
    end do

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, Grid % Comm % nb_tot
      Grid % Comm % bnd_cell_map(c) = int(c-1, SP)
    end do

  !-----------------------!
  !                       !
  !   For parallel runs   !
  !                       !
  !-----------------------!
  else

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!

    !---------------------!
    !   Inside cell map   !
    !---------------------!
    do c = 1, Grid % Comm % nc_sub
      ! Take cell mapping to be the same as global cell numbers but start from 0
      Grid % Comm % cell_map(c) = int(Grid % Comm % cell_glo(c)-1, SP)
    end do

    ! Maps must be in increasing order
    do c = 2, Grid % Comm % nc_sub
      Assert(Grid % Comm % cell_map(c) .gt. Grid % Comm % cell_map(c-1))
    end do

    !-----------------------!
    !   Boundary cell map   !
    !-----------------------!
    cnt = 0
    do reg = Boundary_Regions()
      do c = Cells_In_Region(reg)
        cnt = cnt + 1
        ! Keep in mind there that for boundary cells, cell_glo is negative, so
        ! adding the nb_tot to all the indices, shifts them all to positive
        Grid % Comm % bnd_cell_map(cnt) = int(  Grid % Comm % cell_glo(c)  &
                                              + Grid % Comm % nb_tot, SP)
      end do
    end do  ! region

    ! Maps must be in increasing order
    do c = 2, Grid % Comm % nb_sub
     Assert(Grid % Comm % bnd_cell_map(c) .gt. Grid % Comm % bnd_cell_map(c-1))
    end do

    ! If domain has zero boundary cells, make the only
    ! (fictitious) member in the map point it to zero.
    if(cnt .eq. 0) then
      Grid % Comm % bnd_cell_map(1) = int(0, SP)
      Grid % Comm % nb_f = 0
      Grid % Comm % nb_l = 0
    end if

  end if

  end subroutine
