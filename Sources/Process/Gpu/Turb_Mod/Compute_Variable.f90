!==============================================================================!
  subroutine Compute_Variable(Turb, Grid, phi, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Turb_Type)         :: Turb
  type(Grid_Type),  target :: Grid
  type(Var_Type),   target :: phi
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Val_Type), pointer :: Aval
  type(Sparse_Con_Type), pointer :: Acon
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), visc_eff(:)
  real                           :: urf
  integer                        :: c
!==============================================================================!

  call Profiler % Start('Compute_Variable')

  call Work % Connect_Real_Cell(visc_eff)

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A(MATRIX_ONE)
  val  => Flow % Nat % A(MATRIX_ONE) % val
  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b


  ! Tolerances and under-relaxations are the same for all components
  urf = phi % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(phi % td_scheme .eq. PARABOLIC) then
      phi_oo => phi % oo
      phi_o => phi % o
      !$acc parallel loop independent  &
      !$acc present(  &
      !$acc   grid_region_f_cell,  &
      !$acc   grid_region_l_cell,  &
      !$acc   phi_oo,  &
      !$acc   phi_o   &
      !$acc )
      do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)  ! all present
        phi_oo(c) = phi_o(c)
      end do
      !$acc end parallel
    end if

    phi_o => phi % o
    phi_n => phi % n
    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   phi_o,  &
    !$acc   phi_n   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)  ! all present
      phi_o(c) = phi_n(c)
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------!
  !   Discretize the energy conservation equations   !
  !--------------------------------------------------!
  call Turb % Form_Variable_Matrix(Grid, phi, Flow, Aval, visc_eff,  &
                                   urf, dt=Flow % dt)

  !------------------------------------------------------------!
  !   Insert proper sources (forces) to turbulence equations   !
  !------------------------------------------------------------!

  ! From boundary conditions
  call Turb % Insert_Variable_Bc(Grid, Flow, phi, visc_eff)

  ! Inertial and advection terms
  call Flow % Add_Inertial_Term (Grid, phi, Flow % density)
  call Flow % Add_Advection_Term(Grid, phi, Flow % density)

  ! Insert cross diffusion terms (computers gradients as well)
  call Flow % Add_Cross_Diffusion_Term(Grid, phi, visc_eff)

  !-----------------------------------------!
  !   Insert source term for the variable   !
  !-----------------------------------------!
  if(phi % name .eq. 'VIS') then
    call Turb % Src_Vis_Spalart_Allmaras(Grid, Flow, Acon, Aval)
  end if

  !-----------------------------------------!
  !      Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_Variable_Matrix)   !
  !-----------------------------------------!

  phi_n => phi % n
  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   b,  &
  !$acc   val,  &
  !$acc   dia,  &
  !$acc   phi_n   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * phi_n(c)
  end do
  !$acc end parallel

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Turb_Variable')
  call Flow % Nat % Cg(Acon, Aval, phi % n, b,    &
                       phi % miter, phi % niter,  &
                       phi % tol,   phi % res)
  call Profiler % Stop('CG_for_Turb_Variable')

  do c = Cells_In_Domain_And_Buffers()
    if(phi % n(c) < 0.0) then
      phi % n(c) = phi % o(c)
    end if
  end do

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("turb_variable_solution",              &
                               scalar_name="turb_variable_solution",  &
                               scalar_cell=phi % n)
    call Grid % Save_Debug_Vtu("turb_variable_source",     &
                               inside_name="turb_variable_source",  &
                               inside_cell=b)
# endif

  ! Print info on the screen
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    if(phi % name .eq. 'VIS')  &
      call Info % Iter_Fill_At(3, 1, phi % name,  &
                               phi % res, phi % niter)
  end if

  call Work % Disconnect_Real_Cell(visc_eff)

  call Profiler % Stop('Compute_Variable')

  end subroutine
