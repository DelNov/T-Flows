!==============================================================================!
  subroutine Compute_Pressure(flow, mult, sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt
  integer, intent(in)           :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: a               ! pressure matrix
  type(Matrix_Type), pointer :: m               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: u_f, v_f, w_f, a12, fs
  real                       :: px_f, py_f, pz_f
  character(SL)              :: solver
  real                       :: p_max, p_min, p_nor, p_nor_c, dt
  real :: TMP
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

  call Cpu_Timer_Mod_Start('Compute_Pressure (without solvers)')

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  p      => flow % p
  pp     => flow % pp
  dt     =  flow % dt
  a      => sol % a
  m      => sol % m
  b      => sol % b % val
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, sol, curr_dt, ini)

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

  ! From control file
  call Control_Mod_Normalization_For_Pressure_Solver(p_nor_c)

  ! Calculate pressure magnitude for normalization of pressure solution
  p_max = -HUGE
  p_min = +HUGE
  do c = 1, grid % n_cells
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
  a % val = 0.0

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !    (Convergence is faster with this)    !
  !-----------------------------------------!
  pp % n = 0.0

  !----------------------------------!
  !   Correct fluxes at boundaries   !
  !----------------------------------!
  call Balance_Volume(flow, mult)

  !---------------------------------------!
  !   Compute volume fluxes at internal   !
  !    faces with Rhie and Chow method    !
  !---------------------------------------!
  call Rhie_And_Chow(flow, mult, sol)

  !------------------------------------------!
  !   Update fluxes at boundaries and fill   !
  !   up source term for pressure equation   !
  !------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Internal fluxes fixed with Rhie and Chow method, just update source
    if(c2 > 0) then

      b(c1) = b(c1) - v_flux % n(s)
      b(c2) = b(c2) + v_flux % n(s)

    ! Side is on the boundary
    else

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then

        v_flux % n(s) = ( u % n(c2) * grid % sx(s)     &
                        + v % n(c2) * grid % sy(s)     &
                        + w % n(c2) * grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.   &
              Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) then

        v_flux % n(s) = ( u % n(c2) * grid % sx(s)     &
                        + v % n(c2) * grid % sy(s)     &
                        + w % n(c2) * grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

        a12 = a % fc(s) * grid % vol(c1) / m % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then

        v_flux % n(s) = ( u % n(c1) * grid % sx(s)     &
                        + v % n(c1) * grid % sy(s)     &
                        + w % n(c1) * grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

        a12 = a % fc(s) * grid % vol(c1) / m % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else  ! it is SYMMETRY
        v_flux % n(s) = 0.0
      end if
    end if

  end do

  !-------------------------------------------------------------!
  !   In case of VOF, surface tension and  gravity correction   !
  !-------------------------------------------------------------!
  if(mult % model .eq. VOLUME_OF_FLUID) then
    call Multiphase_Mod_Vof_Pressure_Correction(mult, sol)
  end if

  !-------------------------------------!
  !   Correct fluxes with body forces   !
  !-------------------------------------!
  ! call Field_Mod_Correct_Fluxes_With_Body_Forces(flow, sol)

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Pressure')
  call Solver_Mod_Cg(sol,           &
                     a,             &
                     pp % n,        &
                     b,             &
                     pp % precond,  &
                     pp % mniter,   &      ! max number of iterations
                     pp % eniter,   &      ! executed number of iterations
                     pp % tol,      &
                     pp % res,      &
                     norm = p_nor)         ! number for normalisation

  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Pressure')

  if (flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
  else
    if (flow % i_corr == flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
    end if
  end if

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = 1, grid % n_cells
    p % n(c) =  p % n(c) + pp % urf * pp % n(c)
  end do
  call Grid_Mod_Exchange_Cells_Real(grid, p % n)

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:grid % n_cells))
  p_min  = minval(p % n(1:grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  call Field_Mod_Grad_Pressure_Correction(flow, pp)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(flow, mult, sol, curr_dt, ini)

  call Cpu_Timer_Mod_Stop('Compute_Pressure (without solvers)')

  end subroutine
