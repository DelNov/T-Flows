!==============================================================================!
  subroutine Eddies_Mod_Gather_Bnd_Cells(eddies)
!------------------------------------------------------------------------------!
!   Allocate memory for eddies                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer, allocatable     :: n_bnd_cells_proc(:)
  integer, allocatable     :: s_bnd_cell_proc(:)
  integer                  :: e, c, n, cnt_bnd_cells
!==============================================================================!

  ! Create an alias for grid
  grid => eddies % pnt_grid

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!

  ! Count boundary cells in this processor
  eddies % n_bnd_cells = 0
  do c = -grid % n_bnd_cells, -1
    if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then
      eddies % n_bnd_cells = eddies % n_bnd_cells + 1
    end if
  end do
  print '(a,a,i6)', ' # Number of boundary cells at ',  &
                    trim(eddies % bc_name),             &
                    eddies % n_bnd_cells

  ! Estimate first boundary cell for each processor
  allocate(n_bnd_cells_proc(n_proc));  n_bnd_cells_proc(:) = 0;
  n_bnd_cells_proc(this_proc) = eddies % n_bnd_cells
  call Comm_Mod_Global_Sum_Int_Array(n_proc, n_bnd_cells_proc)

  allocate(s_bnd_cell_proc(n_proc)); s_bnd_cell_proc(:) = 0
  do n = 2, n_proc
    s_bnd_cell_proc(n) = s_bnd_cell_proc(n-1) + n_bnd_cells_proc(n-1)
  end do
  do n = 1, n_proc
    print *, this_proc, n, s_bnd_cell_proc(n)
  end do

  ! Gather coordinates from all processors
  eddies % n_bnd_cells_glo = eddies % n_bnd_cells
  call Comm_Mod_Global_Sum_Int(eddies % n_bnd_cells_glo)
  allocate(eddies % bnd_xc(eddies % n_bnd_cells_glo));  eddies % bnd_xc(:) = 0.0
  allocate(eddies % bnd_yc(eddies % n_bnd_cells_glo));  eddies % bnd_yc(:) = 0.0
  allocate(eddies % bnd_zc(eddies % n_bnd_cells_glo));  eddies % bnd_zc(:) = 0.0

  cnt_bnd_cells = 0
  do c = -grid % n_bnd_cells, -1
    if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then
      cnt_bnd_cells = cnt_bnd_cells + 1
      eddies % bnd_xc(cnt_bnd_cells + s_bnd_cell_proc(this_proc)) = grid % xc(c)
      eddies % bnd_yc(cnt_bnd_cells + s_bnd_cell_proc(this_proc)) = grid % yc(c)
      eddies % bnd_zc(cnt_bnd_cells + s_bnd_cell_proc(this_proc)) = grid % zc(c)
    end if
  end do
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_xc)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_yc)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_zc)

  !-------------------------------!
  !   Place all eddies randomly   !
  !-------------------------------!
  do e = 1, eddies % n_eddies
    call Eddies_Mod_Place_Eddy(eddies, e)
  end do

  end subroutine
