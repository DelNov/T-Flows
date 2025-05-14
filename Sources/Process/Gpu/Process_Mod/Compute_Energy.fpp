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
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), dens_capa(:), cond_eff(:)
  real                           :: urf
  integer                        :: c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Energy')

  call Work % Connect_Real_Cell(dens_capa, cond_eff)

  ! Fill up the dens_capa array
  !$tf-acc loop begin
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()  ! all present
    dens_capa(c) = Flow % density(c) * Flow % capacity(c)
  end do
  !$tf-acc loop end

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % A % dia
  b   => Flow % Nat % b

  ! Tolerances and under-relaxations are the same for all components
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
  call Process % Form_Energy_Matrix(Grid, Flow, Turb, cond_eff,  &
                                    urf, dt=Flow % dt)

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Energy_Bc(Grid, Flow)

  ! Inertial and advection terms
  call Flow % Add_Inertial_Term (Grid, Flow % t, dens_capa)
  call Flow % Add_Advection_Term(Grid, Flow % t, dens_capa)

  ! Insert cross diffusion terms (computers gradients as well)
  call Flow % Add_Cross_Diffusion_Term(Grid, Flow % t, cond_eff)

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
  call Flow % Nat % Cg(Flow % t % n,                        &
                       Flow % t % miter, Flow % t % niter,  &
                       Flow % t % tol,   Flow % t % res)
  call Profiler % Stop('CG_for_Energy')

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("t_0",           &
                               scalar_name="t", &
                               scalar_cell=flow_t_n)
# endif

  call Info % Iter_Fill_At(1, 6, Flow % t % name,  &
                                 Flow % t % res, Flow % t % niter)

  call Work % Disconnect_Real_Cell(dens_capa, cond_eff)

  call Profiler % Stop('Compute_Energy')

  end subroutine
