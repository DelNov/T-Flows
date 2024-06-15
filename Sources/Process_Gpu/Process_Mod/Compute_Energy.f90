!==============================================================================!
  subroutine Compute_Energy(Process, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Val_Type), pointer :: Mval
  type(Sparse_Con_Type), pointer :: Mcon
  real,                  pointer :: val(:)
  integer,               pointer :: dia(:)
  real,                  pointer :: b(:)
  real                           :: tol, fin_res, urf
  integer                        :: nc, n, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Energy')

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  Mcon => Flow % Nat % C
  if(Flow % Nat % use_one_matrix) then
    Mval => Flow % Nat % A(MATRIX_ONE)
    val  => Flow % Nat % A(MATRIX_ONE) % val
  else
    Mval => Flow % Nat % A(MATRIX_T)
    val  => Flow % Nat % A(MATRIX_T) % val
  end if
  dia  => Flow % Nat % C % dia
  b    => Flow % Nat % b
  nc   =  Flow % pnt_grid % n_cells
  fin_res = 0.0

  ! Tolerances and under-relaxations are the same for all components
  tol = Flow % t % tol
  urf = Flow % t % urf

  !---------------------------------------------------!
  !   Update old values (o) and older than old (oo)   !
  !---------------------------------------------------!
  if(Iter % Current() .eq. 1) then

    if(Flow % t % td_scheme .eq. PARABOLIC) then
      !$acc parallel loop independent
      do c = Cells_In_Domain_And_Buffers()
        t_oo(c) = t_o(c)
      end do
      !$acc end parallel
    end if

    !$acc parallel loop independent
    do c = Cells_In_Domain_And_Buffers()
      t_o(c) = t_n(c)
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------!
  !   Discretize the energy conservation equations   !
  !--------------------------------------------------!

  call Process % Form_System_Matrix(Mcon, Mval, Flow, Grid,           &
                                    Flow % density, Flow % capacity,  &
                                    Flow % conductivity,              &
                                    urf, dt = Flow % dt)

  !----------------------------------------------------------!
  !   Insert proper sources (forces) to momentum equations   !
  !----------------------------------------------------------!

  ! From boundary conditions
  call Process % Insert_Energy_Bc(Flow, Grid)

  ! Inertial and advection terms
  call Process % Add_Inertial_Term(Flow % t, Flow, Grid,  &
                                   Flow % density, Flow % capacity)
  call Process % Add_Advection_Term(Flow % t, Flow, Grid,  &
                                    Flow % density, Flow % capacity)

  !---------------------------------------!
  !     Part 2 of the under-relaxation    !
  !   (Part 1 is in Form_System_Matrix)   !
  !---------------------------------------!

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    b(c) = b(c) + val(dia(c)) * (1.0 - urf) * t_n(c)
  end do
  !$acc end parallel

  !------------------------!
  !   Call linear solver   !
  !------------------------!
  call Profiler % Start('CG_for_Energy')
  call Flow % Nat % Cg(Mcon, Mval, t_n, b, nc, n, tol, fin_res)
  call Profiler % Stop('CG_for_Energy')

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("t_0",           &
                               scalar_name="t", &
                               scalar_cell=t_n)
# endif

  call Info % Iter_Fill_At(1, 6, 'T', fin_res, n)

  call Profiler % Stop('Compute_Energy')

  end subroutine
