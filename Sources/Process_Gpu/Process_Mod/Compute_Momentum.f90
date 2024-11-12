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
  real,      contiguous, pointer :: b(:)
  real,      contiguous, pointer :: ui_n(:), ui_o(:), ui_oo(:)
  real                           :: tol, fin_res, urf
  integer                        :: nb, nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Momentum')

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
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  fin_res = 0.0

  if(comp .eq. 1) then
    ui_n  => flow_u_n
    ui_o  => flow_u_o
    ui_oo => flow_u_oo
  else if(comp .eq. 2) then
    ui_n  => flow_v_n
    ui_o  => flow_v_o
    ui_oo => flow_v_oo
  else if(comp .eq. 3) then
    ui_n  => flow_w_n
    ui_o  => flow_w_o
    ui_oo => flow_w_oo
  end if

  ! Tolerances and under-relaxations are the same for all components
  tol = Flow % u % tol
  urf = Flow % u % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(Flow % u % td_scheme .eq. PARABOLIC) then
      !$acc parallel loop independent  &
      !$acc present(  &
      !$acc   grid_region_f_cell,  &
      !$acc   grid_region_l_cell,  &
      !$acc   ui_oo,  &
      !$acc   ui_o   &
      !$acc )
      do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)  ! all present
        ui_oo(c) = ui_o(c)
      end do
      !$acc end parallel
    end if

    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   ui_o,  &
    !$acc   ui_n   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)  ! all present
      ui_o(c) = ui_n(c)
    end do
    !$acc end parallel
  end if

  !----------------------------------------------------!
  !   Discretize the momentum conservation equations   !
  !----------------------------------------------------!

  ! Once is enough, it is the same for all components
  if(comp .eq. 1) then
    call Process % Form_Momentum_Matrix(Grid, Flow, Turb, Acon, Aval,  &
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
    call Flow % Add_Inertial_Term (Grid, Flow % u, flow_density)
    call Flow % Add_Advection_Term(Grid, Flow % u, flow_density)
  else if(comp .eq. 2) then
    call Flow % Add_Inertial_Term (Grid, Flow % v, flow_density)
    call Flow % Add_Advection_Term(Grid, Flow % v, flow_density)
  else if(comp .eq. 3) then
    call Flow % Add_Inertial_Term (Grid, Flow % w, flow_density)
    call Flow % Add_Advection_Term(Grid, Flow % w, flow_density)
  end if

  ! Pressure force
  call Process % Add_Pressure_Term(Grid, Flow, comp=comp)

  !---------------------------------------!
  !    Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_System_Matrix)   !
  !---------------------------------------!

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   b,  &
  !$acc   val,  &
  !$acc   dia,  &
  !$acc   ui_n   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * ui_n(c)
  end do
  !$acc end parallel

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Momentum')
  call Flow % Nat % Cg(Acon, Aval, ui_n, b, nc, n, tol, fin_res)
  call Profiler % Stop('CG_for_Momentum')

  if(comp.eq.1) call Info % Iter_Fill_At(1, 1, 'U', fin_res, n)
  if(comp.eq.2) call Info % Iter_Fill_At(1, 2, 'V', fin_res, n)
  if(comp.eq.3) call Info % Iter_Fill_At(1, 3, 'W', fin_res, n)

  call Profiler % Stop('Compute_Momentum')

  end subroutine
