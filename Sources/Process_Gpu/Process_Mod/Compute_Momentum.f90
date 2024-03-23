!==============================================================================!
  subroutine Compute_Momentum(Process, Flow, Grid, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Val_Type), pointer :: Mval
  type(Sparse_Con_Type), pointer :: Mcon
  real,                  pointer :: val(:)
  integer,               pointer :: dia(:)
  real,                  pointer :: b(:)
  real,                  pointer :: ui_n(:)
  real                           :: tol, fin_res, urf
  integer                        :: nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Momentum')

  Assert(comp .ge. 1)
  Assert(comp .le. 3)

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  Mcon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Mval => Flow % Nat % A(MATRIX_ONE)
    val  => Flow % Nat % A(MATRIX_ONE) % val
  else
    Mval => Flow % Nat % A(MATRIX_UVW)
    val  => Flow % Nat % A(MATRIX_UVW) % val
  end if
  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b
  nc   =  Flow % pnt_grid % n_cells
  fin_res = 0.0

  if(comp .eq. 1) ui_n => u_n
  if(comp .eq. 2) ui_n => v_n
  if(comp .eq. 3) ui_n => w_n

  ! Tolerances and under-relaxations are the same for all components
  tol = Flow % u % tol
  urf = Flow % u % urf

  !----------------------------------------------------!
  !   Discretize the momentum conservation equations   !
  !----------------------------------------------------!

  ! Once is enough, it is the same for all components
  if(comp .eq. 1) then
    call Process % Form_Momentum_Matrix(Flow, Grid, dt=Flow % dt)
  end if

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!
  call Process % Insert_Momentum_Bc(Flow, Grid, comp=comp)
  call Process % Add_Inertial_Term  (Flow, Grid, comp=comp)
  call Process % Add_Advection_Term (Flow, Grid, comp=comp)
  call Process % Add_Pressure_Term  (Flow, Grid, comp=comp)

  !-----------------------------------------!
  !      Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_Momentum_Matrix)   !
  !-----------------------------------------!

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * ui_n(c)
  end do
  !$acc end parallel

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Momentum')
  call Flow % Nat % Cg(Mcon, Mval, ui_n, b, nc, n, tol, fin_res)
  call Profiler % Stop('CG_for_Momentum')

  if(comp.eq.1) call Info % Iter_Fill_At(1, 1, 'U', fin_res, n)
  if(comp.eq.2) call Info % Iter_Fill_At(1, 2, 'V', fin_res, n)
  if(comp.eq.3) call Info % Iter_Fill_At(1, 3, 'W', fin_res, n)

  call Profiler % Stop('Compute_Momentum')

  end subroutine
