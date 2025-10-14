!==============================================================================!
  subroutine Compute_Pressure(Process, Grid, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: val(:)
  real,    contiguous, pointer :: b(:)
  integer, contiguous, pointer :: row(:), col(:)
  real                         :: urf, p_max, p_min, p_nor, p_nor_c
  integer                      :: c, call_type
  integer, save                :: call_count = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Pressure')

  ! Take some aliases
  val => Flow % Nat % A % val
  row => Flow % Nat % A % row
  col => Flow % Nat % A % col
  b   => Flow % Nat % b
  urf =  Flow % pp % urf

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

  ! From control file
  call Control % Normalization_For_Pressure_Solver(p_nor_c)

  ! Calculate pressure magnitude for normalization of pressure solution
  p_max = -HUGE
  p_min = +HUGE
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    p_max = max(p_max, Flow % p % n(c))
    p_min = min(p_min, Flow % p % n(c))
  end do
  !$tf-acc loop end
  call Global % Max_Real(p_max)
  call Global % Min_Real(p_min)
  ! Normalize pressure with the maximum of pressure difference,
  ! value defined in control file and pressure drops.
  p_nor = max( (p_max-p_min), p_nor_c, abs(Flow % bulk % p_drop_x),  &
                                       abs(Flow % bulk % p_drop_y),  &
                                       abs(Flow % bulk % p_drop_z) )

  !-----------------------------------------------!
  !   Discretize the pressure Poisson equations   !
  !-----------------------------------------------!
  call Process % Form_Pressure_Matrix(Flow, Grid)

  !---------------------------------------------------------------!
  !   Insert proper source (volume source) to pressure equation   !
  !---------------------------------------------------------------!
  call Process % Insert_Volume_Source_For_Pressure(Flow, Grid)

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("bp_0",                &
                               inside_name="vol_src", &
                               inside_cell=b)
# endif

  ! Set singularity to the matrix
  if(.not. Flow % has_pressure) then
    call Linalg % Set_Singular(Grid % n_cells, Flow % Nat % A)
  end if

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  if(Flow % pp % solver .eq. 'cg') then
    call Profiler % Start('CG_for_Pressure')
    call Flow % Nat % Cg(Flow % pp % n(1:Grid % n_cells),       &
                         Flow % pp % miter, Flow % pp % niter,  &
                         Flow % pp % tol,   Flow % pp % res,    &
                         norm = p_nor)
    call Profiler % Stop('CG_for_Pressure')
  else if(Flow % pp % solver .eq. 'rs_amg') then
    call Profiler % Start('AMG_for_Pressure')
    call_type = AMG_RUN_ALL_FOUR_STAGES
    if(call_count .gt. 0) call_type = AMG_SOLVE_AND_REPORT
    call Amg % Amg1r5(val, row, col,                    &
                      Flow % pp % n(1:Grid % n_cells),  &
                      b(1:Grid % n_cells),              &
                      Grid % n_cells,                   &
                      Flow % pp % res,                  &
                      call_type)
    Flow % pp % res   = Amg % Final_Residual()
    Flow % pp % niter = Amg % Performed_Cycles()
    call_count = call_count + 1
    call Profiler % Stop('AMG_for_Pressure')
  end if

  call Info % Iter_Fill_At(1, 4, Flow % pp % name,  &
                                 Flow % pp % res, Flow % pp % niter)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    Flow % p % n(c) = Flow % p % n(c) + urf * Flow % pp % n(c)
  end do
  !$tf-acc loop end

  !---------------------------------------------------------------!
  !   Shift the pressure field so that the median value is zero   !
  !---------------------------------------------------------------!
  p_max = -HUGE
  p_min = +HUGE

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    p_max = max(p_max, Flow % p % n(c))
    p_min = min(p_min, Flow % p % n(c))
  end do
  !$tf-acc loop end

  call Global % Max_Real(p_max)
  call Global % Min_Real(p_min)

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    Flow % p % n(c) = Flow % p % n(c) - 0.5 * (p_max + p_min)
  end do
  !$tf-acc loop end

  ! Update buffers for presssure over all processors
  call Grid % Exchange_Cells_Real(Flow % p % n)

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("pp_0",           &
                               scalar_name="pp", &
                               scalar_cell=Flow % pp % n)
    call Grid % Save_Debug_Vtu("p_0",            &
                               scalar_name="p",  &
                               scalar_cell=Flow % p % n)
# endif

  call Profiler % Stop('Compute_Pressure')

  end subroutine
