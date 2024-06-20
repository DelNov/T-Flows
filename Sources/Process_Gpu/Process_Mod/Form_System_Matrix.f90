!==============================================================================!
  subroutine Form_System_Matrix(Process, Acon, Aval, Flow, Grid,  &
                                coef_a, coef_b, diff_coef, urf, dt, save_v_m)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)           :: Process
  type(Sparse_Con_Type), target :: Acon
  type(Sparse_Val_Type), target :: Aval
  type(Field_Type),      target :: Flow
  type(Grid_Type),   intent(in) :: Grid
  real                          :: coef_a(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: coef_b(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: diff_coef(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: urf
  real,    optional, intent(in) :: dt       !! time step
  logical, optional, intent(in) :: save_v_m
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_System_Matrix')

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!

  ! If each varible uses its own matrix
  if(.not. Flow % Nat % use_one_matrix) then
    if(Flow % Nat % A(MATRIX_UVW) % formed) return
  end if

  val => Aval % val
  dia => Acon % dia
  pos => Acon % pos
  fc  => Acon % fc
  nz  =  Acon % nonzeros

  Assert(urf > 0.0)

  !$acc parallel loop independent  &
  !$acc present(val)
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$acc end parallel

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over cells   !
  !--------------------------------------------------!

  !$acc parallel loop independent                                &
  !$acc present(grid_region_f_cell, grid_region_l_cell,          &
  !$acc         grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
  !$acc         val, pos, diff_coef, fc, dia)
  do c1 = Cells_In_Domain_Gpu()  ! all present

    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)

      ! Coefficients inside the domain
      if(c2 .gt. 0) then
        a12 = Face_Value(s, diff_coef(c1), diff_coef(c2)) * fc(s)
        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a12
        end if
        val(dia(c1)) = val(dia(c1)) + a12

      ! Coefficients at the boundaries
      else
        reg = Grid % region % at_cell(c2)
        if(Grid % region % type(reg) .eq. WALL .or.  &
           Grid % region % type(reg) .eq. INFLOW) then
          a12 = diff_coef(c1) * fc(s)
          val(dia(c1)) = val(dia(c1)) + a12
        end if
      end if
    end do
    !$acc end loop

  end do
  !$acc end parallel

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    !$acc parallel loop independent                        &
    !$acc present(grid_region_f_cell, grid_region_l_cell,  &
    !$acc         grid_vol,                                &
    !$acc         val, dia, coef_a, coef_b)
    do c = Cells_In_Domain_Gpu()  ! all present
      val(dia(c)) = val(dia(c)) + coef_a(c) * coef_b(c) * grid_vol(c) / dt
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  if(present(save_v_m)) then
    if(save_v_m) then
      !$acc parallel loop independent                        &
      !$acc present(grid_region_f_cell, grid_region_l_cell,  &
      !$acc         grid_vol,                                &
      !$acc         flow_v_m, val, dia)
      do c = Cells_In_Domain_Gpu()  ! all present
        flow_v_m(c) = grid_vol(c) / val(dia(c))
      end do
      !$acc end parallel

      ! This call is needed, the above loop goes through inside cells only
      call Grid % Exchange_Inside_Cells_Real(flow_v_m)
    end if
  end if

  !-------------------------------------!
  !   Part 1 of the under-relaxation    !
  !   (Part 2 is in Compute_Momentum)   !
  !-------------------------------------!
  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         val, dia)
  do c = Cells_In_Domain_Gpu()  ! all present
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$acc end parallel

  !-------------------------------!
  !   Mark the matrix as formed   !
  !-------------------------------!
  Aval % formed = .true.

# if T_FLOWS_DEBUG == 1
  allocate(temp(Grid % n_cells));  temp(:) = 0.0
  do c = Cells_In_Domain()  ! this is for debugging, don't do it on GPU
    ! or: temp(c) = val(dia(c))
    ! or: temp(c) = Acon % row(c+1) - Acon % row(c)
    temp(c) = flow_v_m(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(temp)
  call Grid % Save_Debug_Vtu("v_m",              &
                             inside_name="v_m",  &
                             inside_cell=temp)
# endif

  call Profiler % Stop('Form_System_Matrix')

  end subroutine
