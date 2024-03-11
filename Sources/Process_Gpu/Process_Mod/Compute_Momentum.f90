!==============================================================================!
  subroutine Compute_Momentum(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Con_Type), pointer :: Mcon
  type(Sparse_Val_Type), pointer :: Mval
  real,                  pointer :: b(:)
  real                           :: tol
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Compute_Momentum')

  Assert(comp .ge. 1)
  Assert(comp .le. 3)

  ! Take some aliases
  Mcon => Flow % Nat % C
  Mval => Flow % Nat % M
  b    => Flow % Nat % b

  ! The tolerances are the same for all components
  tol = Flow % u % tol

  ! Insert proper sources (forces) to momentum equations
  call Process % Insert_Diffusion_Bc(Flow, comp=comp)
  call Process % Add_Inertial_Term  (Flow, comp=comp)
  call Process % Add_Advection_Term (Flow, comp=comp)
  call Process % Add_Pressure_Term  (Flow, comp=comp)

  ! Call linear solver
  call Profiler % Start('CG_for_Momentum')
  if(comp .eq. 1) call Flow % Nat % Cg(Mcon, Mval, u_n, b, grid_n_cells, tol)
  if(comp .eq. 2) call Flow % Nat % Cg(Mcon, Mval, v_n, b, grid_n_cells, tol)
  if(comp .eq. 3) call Flow % Nat % Cg(Mcon, Mval, w_n, b, grid_n_cells, tol)
  call Profiler % Stop('CG_for_Momentum')

  call Profiler % Stop('Compute_Momentum')

  end subroutine
