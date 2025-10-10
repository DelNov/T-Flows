!==============================================================================!
  subroutine Compute_Scalars(Process, Grid, Flow, Turb)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),        pointer :: phi
  real,      contiguous, pointer :: val(:)
  integer,   contiguous, pointer :: dia(:)
  real,      contiguous, pointer :: b(:), diff_eff(:)
  real                           :: urf
  integer                        :: c, sc, row, col
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Scalars')

  call Work % Connect_Real_Cell(diff_eff)

  !-----------------------------!
  !   First take some aliases   !
  !-----------------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % A % dia
  b   => Flow % Nat % b

  !--------------------------------!
  !                                !
  !   Browse through all scalars   !
  !                                !
  !--------------------------------!
  do sc = 1, Flow % n_scalars

    !-------------------------------------------------!
    !   Important: take the alias of current scalar   !
    !-------------------------------------------------!
    phi               => Flow % scalar(sc)
    phi_bnd_cond_type => Flow % scalar(sc) % bnd_cond_type

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
      do c = Cells_In_Domain_And_Buffers() !all present
        phi % o(c) = phi % n(c)
      end do
      !$tf-acc loop end

    end if

    !-------------------------------------------------!
    !   Discretize the scalar conservation equation   !
    !-------------------------------------------------!
    call Process % Form_Scalars_Matrix(Grid, Flow, Turb, diff_eff,  &
                                       sc, urf, dt=Flow % dt)

    !-----------------------------------------------------------!
    !   Insert proper sources (forces) to transport equations   !
    !-----------------------------------------------------------!

    ! From boundary conditions
    call Process % Insert_Scalars_Bc(Grid, Flow, sc)

    ! Inertial and advection terms
    call Flow % Add_Inertial_Term (Grid, phi, Flow % density)
    call Flow % Add_Advection_Term(Grid, phi, Flow % density)

    ! Insert cross diffusion terms (computers gradients as well)
    call Flow % Add_Cross_Diffusion_Term(Grid, phi, diff_eff)

    call User_Mod_Source(Grid, Flow, Turb, phi, sc, Flow % Nat % A, b)

    !----------------------------------------!
    !     Part 2 of the under-relaxation     !
    !   (Part 1 is in Form_Scalars_Matrix)   !
    !----------------------------------------!
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present
      b(c) = b(c) + val(dia(c)) * (1.0 - urf) * phi % n(c)
    end do
    !$tf-acc loop end

    !------------------------!
    !   Call linear solver   !
    !------------------------!
    call Profiler % Start('CG_for_Scalars')
    call Flow % Nat % Cg(phi % n(1:Grid % n_cells),  &
                         phi % miter, phi % niter,   &
                         phi % tol,   phi % res)
    call Profiler % Stop('CG_for_Scalars')

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("C_00",           &
                               scalar_name="C_00", &
                               scalar_cell=phi % n)
# endif

    row = (sc-1)/6 + 1      ! will be 1 (scal. 1-6), 2 (scal. 6-12), etc.
    col = mod(sc-1, 6) + 1  ! will be in range 1 - 6

    call Info % Iter_Fill_Scalar_At(row, col, phi % name,  &
                                              phi % res,   &
                                              phi % niter)
  end do  ! through scalars

  call Work % Disconnect_Real_Cell(diff_eff)

  call Profiler % Stop('Compute_Scalars')

  end subroutine
