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
  val => Flow % Nat % A % val
  dia => Flow % Nat % C % dia
  b   => Flow % Nat % b

  ! Tolerances and under-relaxations are the same for all components
  urf = phi % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(phi % td_scheme .eq. PARABOLIC) then
      !$tf-acc loop begin
      do c = Cells_In_Domain_And_Buffers()  ! all present
        phi % oo(c) = phi % o(c)
      end do
      !$tf-acc loop end
    end if

    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()  ! all present
      phi % o(c) = phi % n(c)
    end do
    !$tf-acc loop end
  end if

  !--------------------------------------------------!
  !   Discretize the energy conservation equations   !
  !--------------------------------------------------!
  call Turb % Form_Variable_Matrix(Grid, phi, Flow, visc_eff,  &
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
    call Turb % Src_Vis_Spalart_Allmaras(Grid, Flow)
  end if

  !-----------------------------------------!
  !      Part 2 of the under-relaxation     !
  !   (Part 1 is in Form_Variable_Matrix)   !
  !-----------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * phi % n(c)
  end do
  !$tf-acc loop end

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Turb_Variable')
  call Flow % Nat % Cg(phi % n,                   &
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
