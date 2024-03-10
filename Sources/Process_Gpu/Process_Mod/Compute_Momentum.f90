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
  integer                    :: n, s, c1, c2, reg
# if T_FLOWS_DEBUG == 1
    character(4) :: name_b = 'bX_Y'
    character(3) :: name_u = 'X_Y'
# endif
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

  ! Check
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      Assert(c2 .lt. 0)
    end do
  end do

  ! Still on aliases
  if(comp .eq. 1) ui_n => Flow % u % n
  if(comp .eq. 2) ui_n => Flow % v % n
  if(comp .eq. 3) ui_n => Flow % w % n
  tol = Flow % u % tol  ! they are the same for all components

  ! Insert proper sources (forces) to momentum equations
  call Process % Insert_Diffusion_Bc(Flow, comp=comp)
  call Process % Add_Inertial_Term  (Flow, comp=comp)
  call Process % Add_Advection_Term (Flow, comp=comp)
  call Process % Add_Pressure_Term  (Flow, comp=comp)

# if T_FLOWS_DEBUG == 1
    write(name_b(4:4), '(i1)') Iter % Current()
    if(comp .eq. 1) write(name_b(2:2), '(a1)') 'u'
    if(comp .eq. 2) write(name_b(2:2), '(a1)') 'v'
    if(comp .eq. 3) write(name_b(2:2), '(a1)') 'w'
    call Grid % Save_Debug_Vtu(name_b, inside_name=name_b, inside_cell=b)
# endif

  ! Call linear solver
  call Profiler % Start('CG_for_Momentum')
  call Flow % Nat % Cg(M, ui_n, b, n, tol)
  call Profiler % Stop('CG_for_Momentum')

# if T_FLOWS_DEBUG == 1
    write(name_u(3:3), '(i1)') Iter % Current()
    if(comp .eq. 1) write(name_u(1:1), '(a1)') 'u'
    if(comp .eq. 2) write(name_u(1:1), '(a1)') 'v'
    if(comp .eq. 3) write(name_u(1:1), '(a1)') 'w'
    call Grid % Save_Debug_Vtu(name_u, scalar_name=name_u, scalar_cell=ui_n)
# endif

  call Profiler % Stop('Compute_Momentum')

  end subroutine
