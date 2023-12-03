!==============================================================================!
  subroutine Gather_Bnd_Cells(Eddies)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to gather boundary cells corresponding to a
!>  specific boundary condition (where eddies are introduced) across all
!>  processors. This subroutine ensures that the data regarding eddy placement
!>  is consistent and continuous across processor boundaries, preventing
!>  discontinuities in the flow field.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias creation: Creates aliases for the grid (Grid) and the flow field   !
!     (Flow) for easier reference.                                             !
!   * Boundary cell counting: Counts the number of boundary cells in the       !
!     current processor that match the specified boundary condition.           !
!     This count is stored in eddies % n_bnd_cells.                            !
!   * Global boundary cell count: If running in a parallel environment, the    !
!     subroutine computes the total number of boundary cells across all        !
!     processors (eddies % n_bnd_cells_glo) and the start index of boundary    !
!     cells for each processor (s_bnd_cell_proc).                              !
!   * Boundary cell data gathering: Collects coordinates (x, y, z), wall       !
!     distance, and flow velocity (u, v, w) for each boundary cell that        !
!     matches the boundary condition. This data is stored in arrays within     !
!     the eddies structure.                                                    !
!   * Inlet plane position: Determines the average position of the inlet       !
!     plane based on the spread (difference between maximum and minimum        !
!     values) of the boundary cell coordinates. This helps in identifying      !
!     which plane (x, y, or z) the eddies are associated with.                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type), target :: Eddies  !! parent class; collection of eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer, allocatable      :: n_bnd_cells_proc(:)
  integer, allocatable      :: s_bnd_cell_proc(:)
  integer                   :: c, n, cnt
  real                      :: x_dff, y_dff, z_dff, x_avg, y_avg, z_avg
!==============================================================================!

  ! Create an alias for Grid
  Grid => Eddies % pnt_grid
  Flow => Eddies % pnt_flow

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!

  ! Count boundary cells in this processor
  Eddies % n_bnd_cells = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      Eddies % n_bnd_cells = Eddies % n_bnd_cells + 1
    end if
  end do

  ! Estimate first boundary cell for each processor
  if(Parallel_Run()) then
    allocate(n_bnd_cells_proc(N_Procs()));  n_bnd_cells_proc(:) = 0;
    n_bnd_cells_proc(This_Proc()) = Eddies % n_bnd_cells
    call Global % Sum_Int_Array(N_Procs(), n_bnd_cells_proc)

    allocate(s_bnd_cell_proc(N_Procs())); s_bnd_cell_proc(:) = 0
    do n = 2, N_Procs()
      s_bnd_cell_proc(n) = s_bnd_cell_proc(n-1) + n_bnd_cells_proc(n-1)
    end do
  else
    allocate(s_bnd_cell_proc(0:0));  s_bnd_cell_proc(:) = 0;
  end if

  ! Gather coordinates from all processors
  Eddies % n_bnd_cells_glo = Eddies % n_bnd_cells
  call Global % Sum_Int(Eddies % n_bnd_cells_glo)
  if(First_Proc()) then
    print '(a,a,i6)', ' # Number of boundary cells at ',  &
                      trim(Eddies % bc_name),             &
                      Eddies % n_bnd_cells_glo
  end if
  allocate(Eddies % bnd_xc(Eddies % n_bnd_cells_glo));  Eddies % bnd_xc(:) = 0.0
  allocate(Eddies % bnd_yc(Eddies % n_bnd_cells_glo));  Eddies % bnd_yc(:) = 0.0
  allocate(Eddies % bnd_zc(Eddies % n_bnd_cells_glo));  Eddies % bnd_zc(:) = 0.0
  allocate(Eddies % bnd_wd(Eddies % n_bnd_cells_glo));  Eddies % bnd_wd(:) = 0.0
  allocate(Eddies % bnd_u (Eddies % n_bnd_cells_glo));  Eddies % bnd_u (:) = 0.0
  allocate(Eddies % bnd_v (Eddies % n_bnd_cells_glo));  Eddies % bnd_v (:) = 0.0
  allocate(Eddies % bnd_w (Eddies % n_bnd_cells_glo));  Eddies % bnd_w (:) = 0.0

  cnt = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      cnt = cnt + 1
      Eddies % bnd_xc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % xc(c)
      Eddies % bnd_yc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % yc(c)
      Eddies % bnd_zc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % zc(c)
      Eddies % bnd_wd(cnt + s_bnd_cell_proc(This_Proc())) = Grid % wall_dist(c)
      Eddies % bnd_u (cnt + s_bnd_cell_proc(This_Proc())) = Flow % u % b(c)
      Eddies % bnd_v (cnt + s_bnd_cell_proc(This_Proc())) = Flow % v % b(c)
      Eddies % bnd_w (cnt + s_bnd_cell_proc(This_Proc())) = Flow % w % b(c)
    end if
  end do
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_xc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_yc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_zc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_wd)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_u)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_v)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_w)

  !-------------------------------------------------------!
  !           Find position of the inlet plane            !
  !   (Important for checking eddies leaving the plane)   !
  !-------------------------------------------------------!
  Eddies % x_plane = HUGE
  Eddies % y_plane = HUGE
  Eddies % z_plane = HUGE

  x_dff = maxval(Eddies % bnd_xc(:)) - minval(Eddies % bnd_xc(:))
  y_dff = maxval(Eddies % bnd_yc(:)) - minval(Eddies % bnd_yc(:))
  z_dff = maxval(Eddies % bnd_zc(:)) - minval(Eddies % bnd_zc(:))

  x_avg = 0.5*(maxval(Eddies % bnd_xc(:)) + minval(Eddies % bnd_xc(:)))
  y_avg = 0.5*(maxval(Eddies % bnd_yc(:)) + minval(Eddies % bnd_yc(:)))
  z_avg = 0.5*(maxval(Eddies % bnd_zc(:)) + minval(Eddies % bnd_zc(:)))

  if(x_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % x_plane = x_avg
  if(y_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % y_plane = y_avg
  if(z_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % z_plane = z_avg

  end subroutine
