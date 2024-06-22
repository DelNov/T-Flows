!==============================================================================!
  subroutine Compute_Energy(Process, Turb, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Turb_Type)          :: Turb
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Val_Type), pointer :: Aval
  type(Sparse_Con_Type), pointer :: Acon
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), dens_capa(:)
  real                           :: tol, fin_res, urf
  integer                        :: nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Energy')

  call Work % Connect_Real_Cell(dens_capa)

  ! Fill up the dens_capa array
  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         dens_capa, flow_density, flow_capacity)
  do c = Cells_At_Boundaries_In_Domain_And_Buffers_Gpu()  ! all present
    dens_capa(c) = flow_density(c) * flow_capacity(c)
  end do
  !$acc end parallel

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  Acon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Aval => Flow % Nat % A(MATRIX_ONE)
    val  => Flow % Nat % A(MATRIX_ONE) % val
  else
    Aval => Flow % Nat % A(MATRIX_T)
    val  => Flow % Nat % A(MATRIX_T) % val
  end if
  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b
  nc   =  Grid % n_cells
  fin_res = 0.0

  ! Tolerances and under-relaxations are the same for all components
  tol = Flow % t % tol
  urf = Flow % t % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(Flow % t % td_scheme .eq. PARABOLIC) then
      !$acc parallel loop independent                        &
      !$acc present(grid_region_f_cell, grid_region_l_cell,  &
      !$acc         flow_t_n, flow_t_o)
      do c = Cells_In_Domain_And_Buffers_Gpu()  ! all present
        flow_t_oo(c) = flow_t_o(c)
      end do
      !$acc end parallel
    end if

    !$acc parallel loop independent                        &
    !$acc present(grid_region_f_cell, grid_region_l_cell,  &
    !$acc         flow_t_n, flow_t_o)
    do c = Cells_In_Domain_And_Buffers_Gpu()  ! all present
      flow_t_o(c) = flow_t_n(c)
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------!
  !   Discretize the energy conservation equations   !
  !--------------------------------------------------!
  call Process % Form_Energy_Matrix(Turb, Flow, Grid, Acon, Aval,  &
                                    urf, dt=Flow % dt)

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Energy_Bc(Flow, Grid)

  ! Inertial and advection terms
  call Process % Add_Inertial_Term (Flow, Grid, Flow % t, dens_capa)
  call Process % Add_Advection_Term(Flow, Grid, Flow % t, dens_capa)

  !---------------------------------------!
  !     Part 2 of the under-relaxation    !
  !   (Part 1 is in Form_Energy_Matrix)   !
  !---------------------------------------!

  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         b, val, dia, flow_t_n)
  do c = Cells_In_Domain_Gpu()  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * flow_t_n(c)
  end do
  !$acc end parallel

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Energy')
  call Flow % Nat % Cg(Acon, Aval, flow_t_n, b, nc, n, tol, fin_res)
  call Profiler % Stop('CG_for_Energy')

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("t_0",           &
                               scalar_name="t", &
                               scalar_cell=flow_t_n)
# endif

  call Info % Iter_Fill_At(1, 6, 'T', fin_res, n)

  call Work % Disconnect_Real_Cell(dens_capa)

  call Profiler % Stop('Compute_Energy')

  end subroutine
