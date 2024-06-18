!==============================================================================!
  subroutine Compute_Pressure(Process, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval
  real, contiguous,      pointer :: b(:)
  real                           :: tol, fin_res, urf, p_max, p_min
  integer                        :: nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Pressure')

  ! Take some aliases
  Acon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Aval => Flow % Nat % A(MATRIX_ONE)
  else
    Aval => Flow % Nat % A(MATRIX_PP)
  end if
  b   => Flow % Nat % b
  tol =  Flow % pp % tol
  urf =  Flow % pp % urf
  nc  =  Grid % n_cells

  !-----------------------------------------------!
  !   Discretize the pressure Poisson equations   !
  !-----------------------------------------------!
  call Process % Form_Pressure_Matrix(Acon, Aval, Flow, Grid)

  !---------------------------------------------------------------!
  !   Insert proper source (volume source) to pressure equation   !
  !---------------------------------------------------------------!
  call Process % Insert_Volume_Source_For_Pressure(Flow, Grid)

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("bp_0",                &
                               inside_name="vol_src", &
                               inside_cell=b)
# endif

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Pressure')
  call Flow % Nat % Cg(Acon, Aval, Flow % pp % n, b, nc, n, tol, fin_res)
  call Profiler % Stop('CG_for_Pressure')

  call Info % Iter_Fill_At(1, 4, 'PP', fin_res, n)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    Flow % p % n(c) = Flow % p % n(c) + urf * Flow % pp % n(c)
  end do
  !$acc end parallel

  !---------------------------------------------------------------!
  !   Shift the pressure field so that the median value is zero   !
  !---------------------------------------------------------------!
  p_max = -HUGE
  p_min = +HUGE

  !$acc parallel loop reduction(max:p_max) reduction(min:p_min)
  do c = Cells_In_Domain()
    p_max = max(p_max, Flow % p % n(c))
    p_min = min(p_min, Flow % p % n(c))
  end do

  call Global % Max_Real(p_max)
  call Global % Min_Real(p_min)

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    Flow % p % n(c) = Flow % p % n(c) - 0.5 * (p_max + p_min)
  end do
  !$acc end parallel

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
