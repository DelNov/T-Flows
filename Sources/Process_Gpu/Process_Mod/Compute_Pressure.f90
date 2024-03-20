!==============================================================================!
  subroutine Compute_Pressure(Proc, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),       pointer :: Grid
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval
  real, contiguous,      pointer :: pp_n(:), b(:)
  real                           :: tol, fin_res
  integer                        :: m, n
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Compute_Pressure')

  ! Take some aliases
  Grid => Flow % pnt_grid
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A
  pp_n => Flow % pp % n
  b    => Flow % Nat % b
  m    =  Flow % pnt_grid % n_cells
  tol  =  Flow % pp % tol

  !---------------------------------------------------------------!
  !   Insert proper source (volume source) to pressure equation   !
  !---------------------------------------------------------------!
  call Process % Insert_Volume_Source_For_Pressure(Flow)

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("bp_0",               &
                               inside_name="vol_src", &
                               inside_cell=b)
# endif

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Pressure')
  call Flow % Nat % Cg(Acon, Aval, pp_n, b, m, n, tol, fin_res)
  call Profiler % Stop('CG_for_Pressure')

  call Info % Iter_Fill_At(1, 4, 'PP', fin_res, n)

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("pp_0",          &
                               scalar_name="pp", &
                               scalar_cell=pp_n)
# endif

  call Profiler % Stop('Compute_Pressure')

  end subroutine
