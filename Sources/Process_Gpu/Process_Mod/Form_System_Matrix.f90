!==============================================================================!
  subroutine Form_System_Matrix(Process, Mcon, Mval, Flow, Grid,  &
                                coef_a, coef_b, diff_coef, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)           :: Process
  type(Sparse_Con_Type), target :: Mcon
  type(Sparse_Val_Type), target :: Mval
  type(Field_Type),      target :: Flow
  type(Grid_Type),   intent(in) :: Grid
  real                          :: coef_a(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: coef_b(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: diff_coef(-Grid % n_bnd_cells:Grid % n_cells)
  real                          :: urf
  real,  optional,   intent(in) :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), v_m(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12
# if T_FLOWS_DEBUG == 1
  real, allocatable :: work(:)
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

  val => Mval % val
  dia => Mcon % dia
  pos => Mcon % pos
  fc  => Mcon % fc
  nz  =  Mcon % nonzeros
  v_m => Flow % v_m

  Assert(urf > 0.0)

  !$acc kernels
  do i = 1, nz
    val(i) = 0.0
  end do
  !$acc end kernels

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over cells   !
  !--------------------------------------------------!

  !$acc parallel loop independent
  do c1 = Cells_In_Domain()  ! that should be sufficient

    !$acc loop seq
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      ! Coefficients inside the domain
      if(c2 .gt. 0) then
        a12 = 0.5 * (diff_coef(c1)+diff_coef(c2)) * fc(s)
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
    !$acc parallel loop independent
    do c = Cells_In_Domain()
      val(dia(c)) = val(dia(c)) + coef_a(c) * coef_b(c) * Grid % vol(c) / dt
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    v_m(c) = Grid % vol(c) / val(dia(c))
  end do
  !$acc end parallel

  call Grid % Exchange_Inside_Cells_Real(v_m)

  !-------------------------------------!
  !   Part 1 of the under-relaxation    !
  !   (Part 2 is in Compute_Momentum)   !
  !-------------------------------------!
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$acc end parallel

  !-------------------------------!
  !   Mark the matrix as formed   !
  !-------------------------------!
  Mval % formed = .true.

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = Cells_In_Domain()
    ! or: work(c) = val(dia(c))
    ! or: work(c) = Mcon % row(c+1) - Mcon % row(c)
    work(c) = v_m(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(work)
  call Grid % Save_Debug_Vtu("v_m",              &
                             inside_name="v_m",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_System_Matrix')

  end subroutine
