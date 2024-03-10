!==============================================================================!
  subroutine Compute_Momentum(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Sparse_Type), pointer :: M
  real,              pointer :: ui_n(:)
  real,              pointer :: b(:)
  real                       :: tol
  integer                    :: n
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Compute_Momentum')

  Assert(comp .ge. 1)
  Assert(comp .le. 3)

  ! Take some aliases
  Grid => Flow % pnt_grid
  M    => Flow % Nat % M
  b    => Flow % Nat % b
  n    =  Grid % n_cells

  ! Still on aliases
  if(comp .eq. 1) ui_n => Flow % u % n
  if(comp .eq. 2) ui_n => Flow % v % n
  if(comp .eq. 3) ui_n => Flow % w % n
  tol = Flow % u % tol  ! they are the same for all components

  ! Insert proper sources (forces) to momentum equations
  call Process % Insert_Diffusion_Bc(Flow, comp=comp)
  call Process % Add_Inertial_Term  (Flow, comp=comp)
  call Process % Add_Advection_Term (Flow, comp=comp)
  call Process % Add_Body_Force_Term(Flow, comp=comp)
  call Process % Add_Pressure_Term  (Flow, comp=comp)

# if T_FLOWS_DEBUG == 1
  !@ if(comp .eq. 1) call Grid % Save_Debug_Vtu("bu_0",                 &
  !@                                             inside_name="u_force", &
  !@                                             inside_cell=b)
  !@ if(comp .eq. 2) call Grid % Save_Debug_Vtu("bv_0",                 &
  !@                                             inside_name="v_force", &
  !@                                             inside_cell=b)
  !@ if(comp .eq. 3) call Grid % Save_Debug_Vtu("bw_0",                 &
  !@                                             inside_name="w_force", &
  !@                                             inside_cell=b)
# endif

  ! Call linear solver
  call Profiler % Start('CG_for_Momentum')
  call Flow % Nat % Cg(M, ui_n, b, n, tol)
  call Profiler % Stop('CG_for_Momentum')

# if T_FLOWS_DEBUG == 1
  !@ if(comp .eq. 1) call Grid % Save_Debug_Vtu("u_0",                 &
  !@                                             scalar_name="u_comp", &
  !@                                             scalar_cell=ui_n)
  !@ if(comp .eq. 2) call Grid % Save_Debug_Vtu("v_0",                 &
  !@                                             scalar_name="v_comp", &
  !@                                             scalar_cell=ui_n)
  !@ if(comp .eq. 3) call Grid % Save_Debug_Vtu("w_0",                 &
  !@                                             scalar_name="w_comp", &
  !@                                             scalar_cell=ui_n)
# endif

  call Profiler % Stop('Compute_Momentum')

  end subroutine
