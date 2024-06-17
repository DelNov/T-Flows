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
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), ones(:)
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
  Mcon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Mval => Flow % Nat % A(MATRIX_ONE)
    val  => Flow % Nat % A(MATRIX_ONE) % val
  else
    Mval => Flow % Nat % A(MATRIX_UVW)
    val  => Flow % Nat % A(MATRIX_UVW) % val
  end if

  ones => Flow % work
  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b
  nb   =  Grid % n_bnd_cells
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

  ! Set work variable (Buoyancy_Forces uses it!)

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    ones(c) = 1.0
  end do
  !$acc end parallel

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(Flow % u % td_scheme .eq. PARABOLIC) then
      !$acc parallel loop independent
      do c = Cells_In_Domain_And_Buffers()
        ui_oo(c) = ui_o(c)
      end do
      !$acc end parallel
    end if

    !$acc parallel loop independent
    do c = Cells_In_Domain_And_Buffers()
      ui_o(c) = ui_n(c)
    end do
    !$acc end parallel
  end if

  !----------------------------------------------------!
  !   Discretize the momentum conservation equations   !
  !----------------------------------------------------!

  ! Once is enough, it is the same for all components
  if(comp .eq. 1) then
    call Process % Form_System_Matrix(Mcon, Mval, Flow, Grid,  &
                                      Flow % density, ones,    &
                                      Flow % viscosity,        &
                                      urf, dt = Flow % dt)
  end if

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Momentum_Bc(Flow, Grid, comp=comp)

  ! Buoyancy forces
  call Flow % Buoyancy_Forces(Grid, comp)

  ! Set work variable (Buoyancy_Forces uses it!)
  !$acc parallel loop independent
  do c = Cells_In_Domain_And_Buffers()
    ones(c) = 1.0
  end do
  !$acc end parallel

  ! Inertial and advection terms
  if(comp .eq. 1) then
    call Process % Add_Inertial_Term(Flow % u, Flow, Grid,  &
                                     Flow % density, ones)
    call Process % Add_Advection_Term(Flow % u, Flow, Grid,  &
                                      Flow % density, ones)
  else if(comp .eq. 2) then
    call Process % Add_Inertial_Term(Flow % v, Flow, Grid,  &
                                     Flow % density, ones)
    call Process % Add_Advection_Term(Flow % v, Flow, Grid,  &
                                      Flow % density, ones)
  else if(comp .eq. 3) then
    call Process % Add_Inertial_Term(Flow % w, Flow, Grid,  &
                                     Flow % density, ones)
    call Process % Add_Advection_Term(Flow % w, Flow, Grid,  &
                                      Flow % density, ones)
  end if

  ! Pressure force
  call Process % Add_Pressure_Term(Flow, Grid, comp=comp)

  !---------------------------------------!
  !    Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_System_Matrix)   !
  !---------------------------------------!

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
