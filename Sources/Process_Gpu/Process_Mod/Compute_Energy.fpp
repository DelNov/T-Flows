!==============================================================================!
  subroutine Compute_Energy(Process, Grid, Flow, Turb)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
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

      !$tf-acc loop begin
      do c = Cells_In_Domain_And_Buffers()  ! all present
        Flow % t % oo(c) = Flow % t % o(c)
      end do
      !$tf-acc loop end

    end if

    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()  ! all present
      Flow % t % o(c) = Flow % t % n(c)
    end do
    !$tf-acc loop end

  end if

  !--------------------------------------------------!
  !   Discretize the energy conservation equations   !
  !--------------------------------------------------!
  call Process % Form_Energy_Matrix(Grid, Flow, Turb, Acon, Aval,  &
                                    urf, dt=Flow % dt)

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Energy_Bc(Grid, Flow)

  ! Inertial and advection terms
  call Process % Add_Inertial_Term (Grid, Flow, Flow % t, dens_capa)
  call Process % Add_Advection_Term(Grid, Flow, Flow % t, dens_capa)

  !---------------------------------------!
  !     Part 2 of the under-relaxation    !
  !   (Part 1 is in Form_Energy_Matrix)   !
  !---------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * Flow % t % n(c)
  end do
  !$tf-acc loop end

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
