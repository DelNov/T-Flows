!==============================================================================!
  subroutine Compute_Momentum(Process, Grid, Flow, Turb, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Val_Type), pointer :: Aval
  type(Sparse_Con_Type), pointer :: Acon
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), visc_eff(:)
  real,      contiguous, pointer :: ui_n(:), ui_o(:), ui_oo(:)
  real                           :: tol, fin_res, urf
  integer                        :: nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Momentum')

  call Work % Connect_Real_Cell(visc_eff)

  Assert(comp .ge. 1)
  Assert(comp .le. 3)

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  Acon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Aval => Flow % Nat % A(MATRIX_ONE)
    val  => Flow % Nat % A(MATRIX_ONE) % val
  else
    Aval => Flow % Nat % A(MATRIX_UVW)
    val  => Flow % Nat % A(MATRIX_UVW) % val
  end if

  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b
  nc   =  Grid % n_cells
  fin_res = 0.0

  if(comp .eq. 1) then
    ui_n  => Flow % u % n
    ui_o  => Flow % u % o
    ui_oo => Flow % u % oo
  else if(comp .eq. 2) then
    ui_n  => Flow % v % n
    ui_o  => Flow % v % o
    ui_oo => Flow % v % oo
  else if(comp .eq. 3) then
    ui_n  => Flow % w % n
    ui_o  => Flow % w % o
    ui_oo => Flow % w % oo
  end if

  ! Tolerances and under-relaxations are the same for all components
  tol = Flow % u % tol
  urf = Flow % u % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(Flow % u % td_scheme .eq. PARABOLIC) then
      !$tf-acc loop begin
      do c = Cells_In_Domain_And_Buffers()  ! all present
        ui_oo(c) = ui_o(c)
      end do
      !$tf-acc loop end
    end if

    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()  ! all present
      ui_o(c) = ui_n(c)
    end do
    !$tf-acc loop end
  end if

  !----------------------------------------------------!
  !   Discretize the momentum conservation equations   !
  !----------------------------------------------------!

  ! Once is enough, it is the same for all components
  if(comp .eq. 1) then
    call Process % Form_Momentum_Matrix(Grid, Flow, Turb, Aval, visc_eff,  &
                                        urf, dt = Flow % dt)
  end if

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Momentum_Bc(Grid, Flow, comp=comp)

  ! Buoyancy forces
  call Flow % Buoyancy_Forces(Grid, comp)

  ! Inertial and advection terms
  if(comp .eq. 1) then
    call Flow % Add_Inertial_Term (Grid, Flow % u, Flow % density)
    call Flow % Add_Advection_Term(Grid, Flow % u, Flow % density)
    call Flow % Add_Cross_Diffusion_Term(Grid, Flow % u, visc_eff)
  else if(comp .eq. 2) then
    call Flow % Add_Inertial_Term (Grid, Flow % v, Flow % density)
    call Flow % Add_Advection_Term(Grid, Flow % v, Flow % density)
    call Flow % Add_Cross_Diffusion_Term(Grid, Flow % v, visc_eff)
  else if(comp .eq. 3) then
    call Flow % Add_Inertial_Term (Grid, Flow % w, Flow % density)
    call Flow % Add_Advection_Term(Grid, Flow % w, Flow % density)
    call Flow % Add_Cross_Diffusion_Term(Grid, Flow % w, visc_eff)
  end if

  ! Pressure force
  call Flow % Grad_Pressure(Grid, Flow % p)
  call Process % Add_Pressure_Term(Grid, Flow, comp=comp)

  !---------------------------------------!
  !    Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_System_Matrix)   !
  !---------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * ui_n(c)
  end do
  !$tf-acc loop end

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Momentum')
  call Flow % Nat % Cg(Acon, Aval, ui_n, b,  &
                       Flow % u % miter, n, tol, fin_res)
  call Profiler % Stop('CG_for_Momentum')

  if(comp.eq.1) call Info % Iter_Fill_At(1, 1, Flow % u % name, fin_res, n)
  if(comp.eq.2) call Info % Iter_Fill_At(1, 2, Flow % v % name, fin_res, n)
  if(comp.eq.3) call Info % Iter_Fill_At(1, 3, Flow % w % name, fin_res, n)

  call Work % Disconnect_Real_Cell(visc_eff)

  call Profiler % Stop('Compute_Momentum')

  end subroutine
