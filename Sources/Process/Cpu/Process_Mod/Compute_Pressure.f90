!==============================================================================!
  subroutine Compute_Pressure(Process, Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!>  The subroutine Compute_Pressure in T-Flows is primarily focused on forming
!>  and solving the pressure equation for the SIMPLE method.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and setup: Commences with initializing essential          !
!     variables and pointers, such as those for the grid, flow field, volume   !
!     fluxes, and pressure variables. Also starts the performance profiler.    !
!   * Equation discretization: This subroutine discretizes the pressure        !
!     equation as per the SIMPLE algorithm. It takes into account the volume   !
!     fluxes at boundaries and implements the Rhie and Chow interpolation to   !
!     avoid checkerboard pressure patterns.                                    !
!   * Spatial discretization and matrix assembly: It iterates through all the  !
!     faces of the computational grid to compute coefficients of the           !
!     discretized pressure equation. It assembles the system matrix,           !
!     considering internal and boundary fluxes, and additional sources due to  !
!     mass transfer in VOF simulations.                                        !
!   * Solving equations: The subroutine employs a linear solver to solve the   !
!     discretized pressure correction equation. It normalizes the pressure     !
!     solution to aid convergence and handles the singularity of the pressure  !
!     matrix in incompressible flow simulations.                               !
!   * Post-processing: After solving the equations, updates the pressure field !
!     by adding pressure corrections. It normalizes the updated pressure field !
!     to remove any arbitrary constants, and recalculates the pressure         !
!     gradients for the updated pressure field.                                !
!   * User-defined functions: Invokes various user-defined functions at        !
!     specific stages for additional customization and integration of specific !
!     behaviors into the pressure computation process.                         !
!   * Performance monitoring: Consistently monitors performance to assist in   !
!     the analysis and optimization of the simulation, helping identify        !
!     potential bottlenecks.                                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Vof_Type),      target :: Vof      !! VOF object
  type(Solver_Type),   target :: Sol      !! solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: p_max, p_min, dt, a12
!------------------------------------------------------------------------------!
!
!   The form of equations which are being solved:
!
!      /           /
!     |           |
!     | u dS = dt | GRAD pp dS
!     |           |
!    /           /
!
!   Dimension of the system under consideration
!
!     [A] {pp} = {b}     [m^3/s]
!
!   Dimensions of certain variables
!
!     A                     [m^4s/kg]
!     pp                    [kg/(ms^2)]
!     p % x, p % y, p % z   [kg/(m^2 s^2)]
!     px_f, py_f, pz_f      [kg/(m^2 s^2)]
!     b                     [m^3/s]
!     v_flux                [m^3/s]
!
!==============================================================================!

  call Profiler % Start('Compute_Pressure (without solvers)')

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

  ! Volume balance reporting
  call Flow % Report_Vol_Balance_Start(Iter % Current())

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(Flow, Vof, Sol)

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

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
  call Process % Balance_Volume(Flow, Vof)

  !---------------------------------------!
  !   Compute volume fluxes at internal   !
  !    faces with Rhie and Chow method    !
  !---------------------------------------!
  call Process % Rhie_And_Chow(Flow, Vof, Sol % Nat)

  !-----------------------------------------------!
  !   Update pressure r.h.s. with volume fluxes   !
  !-----------------------------------------------!
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
        a12 = A % fc(s) * M % v_m(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
      end if
    end if
  end do

  ! Volume balance reporting
  call Flow % Report_Vol_Balance(Sol, Iter % Current())

  !------------------------------------------!
  !   Cross diffusion fluxes for pressure    !
  !- - - - - - - - - - - - - - - - - - - - - +-------------------------!
  !   They are correct in the present form, but they have very small   !
  !   impact on results while being quite detremental on convergence   !
  !--------------------------------------------------------------------!
  !@  do s = 1, Grid % n_faces
  !@    c1 = Grid % faces_c(1,s)
  !@    c2 = Grid % faces_c(2,s)
  !@    fs = Grid % f(s)
  !@
  !@    if(c2 > 0) then
  !@
  !@      ! Unit: (m^3 s)/kg
  !@      v_m_c1 = M % v_m(c1)
  !@      v_m_c2 = M % v_m(c2)
  !@      v_m_f  = fs * v_m_c1 + (1.0-fs) * v_m_c2
  !@
  !@      ! Interpolate pressure gradients
  !@      ! Unit: kg/(m^2 s^2)
  !@      px_f = fs * pp % x(c1) + (1.0-fs) * pp % x(c2)
  !@      py_f = fs * pp % y(c1) + (1.0-fs) * pp % y(c2)
  !@      pz_f = fs * pp % z(c1) + (1.0-fs) * pp % z(c2)
  !@
  !@      ! Explicit and implicit pressure correction "fluxes"
  !@      ! Unit: (m^3 s)/kg * kg/(m^2 s^2) * m^2 = m^3/s
  !@      f_ex = v_m_f * (   px_f * Grid % sx(s)   &
  !@                       + py_f * Grid % sy(s)   &
  !@                       + pz_f * Grid % sz(s))
  !@
  !@      ! Unit: (m^3 s)/kg * m * kg/(m^2 s^2) * m = m^3/s
  !@      f_im = v_m_f * A % fc(s) * (   px_f * Grid % dx(s)    &
  !@                                   + py_f * Grid % dy(s)    &
  !@                                   + pz_f * Grid % dz(s) )
  !@    end if
  !@  end do

  !-------------------------------------------------------------------------!
  !   In case of mass transfer, add addtional source to pressure equation   !
  !-------------------------------------------------------------------------!
  call Vof % Mass_Transfer_Pressure_Source(b)

  call Profiler % Start(String % First_Upper(pp % solver)  //  &
                        ' (solver for pressure)')

  ! Set singularity to the matrix
  if(.not. Flow % has_pressure) then
    call Sol % Set_Singular(pp)
  end if

  ! Call linear solver
  call Sol % Run(A, pp, b)

  ! Remove singularity from the matrix
  if(.not. Flow % has_pressure) then
    call Sol % Remove_Singular(pp)
  end if

  call Profiler % Stop(String % First_Upper(pp % solver)  //  &
                       ' (solver for pressure)')

  if (Flow % p_m_coupling == SIMPLE) then
    call Info % Iter_Fill_At(1, 4, pp % name, pp % res, pp % niter)
  else
    if (Flow % i_corr == Flow % n_piso_corrections) then
      call Info % Iter_Fill_At(1, 4, pp % name, pp % res, pp % niter)
    end if
  end if

  call Flow % Grad_Pressure(pp)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()
    p % n(c) =  p % n(c) + pp % urf * pp % n(c)
  end do

  call Flow % Grad_Pressure(p)

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:Grid % n_cells))
  p_min  = minval(p % n(1:Grid % n_cells))

  call Global % Max_Real(p_max)
  call Global % Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol)

  ! Volume balance reporting
  call Flow % Report_Vol_Balance_Stop()

  call Profiler % Stop('Compute_Pressure (without solvers)')

  end subroutine
