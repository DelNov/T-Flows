!==============================================================================!
  subroutine Form_Energy_Matrix(Process, Acon, Aval, Flow, Grid,  &
                                dens_capa, cond, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)           :: Process
  type(Sparse_Con_Type), target :: Acon
  type(Sparse_Val_Type), target :: Aval
  type(Field_Type),      target :: Flow
  type(Grid_Type),   intent(in) :: Grid
  real                          :: dens_capa(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: cond(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: urf
  real,    optional, intent(in) :: dt       !! time step
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

  call Profiler % Start('Form_Energy_Matrix')

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
  !$acc         val, pos, cond, fc, dia)
  do c1 = Cells_In_Domain_Gpu()  ! all present

    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)

      ! Coefficients inside the domain
      if(c2 .gt. 0) then
        a12 = Face_Value(s, cond(c1), cond(c2)) * fc(s)
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
          a12 = cond(c1) * fc(s)
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
    !$acc         val, dia, dens_capa)
    do c = Cells_In_Domain_Gpu()  ! all present
      val(dia(c)) = val(dia(c)) + dens_capa(c) * grid_vol(c) / dt
    end do
    !$acc end parallel
  end if

  !------------------------------------!
  !   Part 1 of the under-relaxation   !
  !   (Part 2 is in Compute_Energy)    !
  !------------------------------------!
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

  call Profiler % Stop('Form_Energy_Matrix')

  end subroutine
