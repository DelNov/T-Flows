!==============================================================================!
  subroutine Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: BEGIN = 12
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  integer, save              :: total_cells
  real                       :: total_source
  real                       :: p_max, p_min, p_nor, p_nor_c, dt, a12
  character(SL)              :: solver
!==============================================================================!
!
!   The form of equations which I am solving:
!
!      /           /
!     |           |
!     | u dS = dt | GRAD pp dS
!     |           |
!    /           /
!
!   Dimension of the system under consideration
!
!     [App] {pp} = {bpp}     [m^3/s]
!
!   Dimensions of certain variables
!
!     app                   [m^4s/kg]
!     pp                    [kg/ms^2]
!     p % x, p % y, p % z   [kg/(m^2 s^2)]
!     px_f, py_f, pz_f      [kg/(m^2 s^2)]
!     b                     [m^3/s]
!     v_flux                [m^3/s]
!
!==============================================================================!

  call Cpu_Timer % Start('Compute_Pressure (without solvers)')

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  p      => Flow % p
  pp     => Flow % pp
  dt     =  Flow % dt
  A      => Sol % Nat % A
  M      => Sol % Nat % M
  b      => Sol % Nat % b % val
  call Flow % Alias_Momentum(u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

  ! From control file
  call Control_Mod_Normalization_For_Pressure_Solver(p_nor_c)

  ! Calculate pressure magnitude for normalization of pressure solution
  p_max = -HUGE
  p_min = +HUGE
  do c = 1, Grid % n_cells
    p_max = max(p_max, p % n(c))
    p_min = min(p_min, p % n(c))
  end do
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  ! Normalize pressure with the maximum of pressure difference, 
  ! value defined in control file and pressure drops.
  p_nor = max( (p_max-p_min), p_nor_c, abs(bulk % p_drop_x),  &
                                       abs(bulk % p_drop_y),  &
                                       abs(bulk % p_drop_z) )

  ! Initialize matrix and right hand side
  b       = 0.0
  A % val = 0.0

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !    (Convergence is faster with this)    !
  !-----------------------------------------!
  pp % n = 0.0

  !----------------------------------!
  !   Correct fluxes at boundaries   !
  !----------------------------------!
  call Balance_Volume(Flow, Vof)

  !---------------------------------------!
  !   Compute volume fluxes at internal   !
  !    faces with Rhie and Chow method    !
  !---------------------------------------!
  call Rhie_And_Chow(Flow, Vof, Sol)

  b(:) = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Internal fluxes fixed with Rhie and Chow method, just update source
    if(c2 > 0) then

      b(c1) = b(c1) - v_flux % n(s)
      b(c2) = b(c2) + v_flux % n(s)

    ! Side is on the boundary
    else

      b(c1) = b(c1) - v_flux % n(s)

      if(Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then
        a12 = A % fc(s) * Grid % vol(c1) / M % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
      end if
    end if

  end do

  !-------------------------------------------------------------------------!
  !   In case of mass transfer, add addtional source to pressure equation   !
  !-------------------------------------------------------------------------!
  call Vof % Mass_Transfer_Pressure_Source(b)

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer % Start('Linear_Solver_For_Pressure')

  ! Set singularity to the matrix
  if(.not. Flow % has_pressure) then
    call Sol % Set_Singular(A)
  end if

  ! Call linear solver
  call Sol % Run(pp % solver,     &
                 pp % prec,       &
                 pp % prec_opts,  &
                 A,               &
                 pp % n,          &
                 b,               &
                 pp % mniter,     &  ! max number of iterations
                 pp % eniter,     &  ! executed number of iterations
                 pp % tol,        &  ! tolerance
                 pp % res,        &  ! final residual
                 norm = p_nor)       ! number for normalisation

  ! Remove singularity from the matrix
  if(.not. Flow % has_pressure) then
    call Sol % Remove_Singular(A)
  end if

  call Cpu_Timer % Stop('Linear_Solver_For_Pressure')

  if (Flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
  else
    if (Flow % i_corr == Flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
    end if
  end if

  call Flow % Grad_Pressure_Correction(pp)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    p % n(c) =  p % n(c) + pp % urf * pp % n(c)
  end do

  call Flow % Grad_Pressure(p)

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:Grid % n_cells))
  p_min  = minval(p % n(1:Grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)

  call Cpu_Timer % Stop('Compute_Pressure (without solvers)')

  end subroutine
