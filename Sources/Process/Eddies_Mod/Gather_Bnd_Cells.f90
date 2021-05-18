!==============================================================================!
  subroutine Eddies_Mod_Gather_Bnd_Cells(eddies)
!------------------------------------------------------------------------------!
!   Gather all boundary cells where eddies are introduced                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer, allocatable      :: n_bnd_cells_proc(:)
  integer, allocatable      :: s_bnd_cell_proc(:)
  integer                   :: c, n, cnt
  real                      :: x_dff, y_dff, z_dff, x_avg, y_avg, z_avg
!==============================================================================!

  ! Create an alias for Grid
  Grid => eddies % pnt_grid
  Flow => eddies % pnt_flow

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!

  ! Count boundary cells in this processor
  eddies % n_bnd_cells = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name( c) .eq. eddies % bc_name) then
      eddies % n_bnd_cells = eddies % n_bnd_cells + 1
    end if
  end do

  ! Estimate first boundary cell for each processor
  if(n_proc > 1) then
    allocate(n_bnd_cells_proc(n_proc));  n_bnd_cells_proc(:) = 0;
    n_bnd_cells_proc(this_proc) = eddies % n_bnd_cells
    call Comm_Mod_Global_Sum_Int_Array(n_proc, n_bnd_cells_proc)

    allocate(s_bnd_cell_proc(n_proc)); s_bnd_cell_proc(:) = 0
    do n = 2, n_proc
      s_bnd_cell_proc(n) = s_bnd_cell_proc(n-1) + n_bnd_cells_proc(n-1)
    end do
  else
    allocate(s_bnd_cell_proc(0:0));  s_bnd_cell_proc(:) = 0;
  end if

  ! Gather coordinates from all processors
  eddies % n_bnd_cells_glo = eddies % n_bnd_cells
  call Comm_Mod_Global_Sum_Int(eddies % n_bnd_cells_glo)
  if(this_proc < 2) then
    print '(a,a,i6)', ' # Number of boundary cells at ',  &
                      trim(eddies % bc_name),             &
                      eddies % n_bnd_cells_glo
  end if
  allocate(eddies % bnd_xc(eddies % n_bnd_cells_glo));  eddies % bnd_xc(:) = 0.0
  allocate(eddies % bnd_yc(eddies % n_bnd_cells_glo));  eddies % bnd_yc(:) = 0.0
  allocate(eddies % bnd_zc(eddies % n_bnd_cells_glo));  eddies % bnd_zc(:) = 0.0
  allocate(eddies % bnd_wd(eddies % n_bnd_cells_glo));  eddies % bnd_wd(:) = 0.0
  allocate(eddies % bnd_u (eddies % n_bnd_cells_glo));  eddies % bnd_u (:) = 0.0
  allocate(eddies % bnd_v (eddies % n_bnd_cells_glo));  eddies % bnd_v (:) = 0.0
  allocate(eddies % bnd_w (eddies % n_bnd_cells_glo));  eddies % bnd_w (:) = 0.0

  cnt = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name( c) .eq. eddies % bc_name) then
      cnt = cnt + 1
      eddies % bnd_xc(cnt + s_bnd_cell_proc(this_proc)) = Grid % xc(c)
      eddies % bnd_yc(cnt + s_bnd_cell_proc(this_proc)) = Grid % yc(c)
      eddies % bnd_zc(cnt + s_bnd_cell_proc(this_proc)) = Grid % zc(c)
      eddies % bnd_wd(cnt + s_bnd_cell_proc(this_proc)) = Grid % wall_dist(c)
      eddies % bnd_u (cnt + s_bnd_cell_proc(this_proc)) = Flow % u % b(c)
      eddies % bnd_v (cnt + s_bnd_cell_proc(this_proc)) = Flow % v % b(c)
      eddies % bnd_w (cnt + s_bnd_cell_proc(this_proc)) = Flow % w % b(c)
    end if
  end do
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_xc)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_yc)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_zc)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_wd)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_u)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_v)
  call Comm_Mod_Global_Sum_Real_Array(eddies % n_bnd_cells_glo, eddies % bnd_w)

  !-------------------------------------------------------!
  !           Find position of the inlet plane            !
  !   (Important for checking eddies leaving the plane)   !
  !-------------------------------------------------------!
  eddies % x_plane = HUGE
  eddies % y_plane = HUGE
  eddies % z_plane = HUGE

  x_dff = maxval(eddies % bnd_xc(:)) - minval(eddies % bnd_xc(:))
  y_dff = maxval(eddies % bnd_yc(:)) - minval(eddies % bnd_yc(:))
  z_dff = maxval(eddies % bnd_zc(:)) - minval(eddies % bnd_zc(:))

  x_avg = 0.5*(maxval(eddies % bnd_xc(:)) + minval(eddies % bnd_xc(:)))
  y_avg = 0.5*(maxval(eddies % bnd_yc(:)) + minval(eddies % bnd_yc(:)))
  z_avg = 0.5*(maxval(eddies % bnd_zc(:)) + minval(eddies % bnd_zc(:)))

  if(x_dff .eq. min(x_dff, y_dff, z_dff)) eddies % x_plane = x_avg
  if(y_dff .eq. min(x_dff, y_dff, z_dff)) eddies % y_plane = y_avg
  if(z_dff .eq. min(x_dff, y_dff, z_dff)) eddies % z_plane = z_avg

  end subroutine
