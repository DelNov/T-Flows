!==============================================================================!
  subroutine Compute_Pressure(flow, mult, sol, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Const_Mod
  use Cpu_Timer_Mod,  only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use Field_Mod
  use Grid_Mod,       only: Grid_Type
  use Bulk_Mod,       only: Bulk_Type
  use Info_Mod,       only: Info_Mod_Iter_Fill_At
  use Solver_Mod,     only: Solver_Type, Bicg, Cg, Cgs, Acm
  use Matrix_Mod,     only: Matrix_Type
  use Control_Mod
  use Multiphase_Mod, only: Multiphase_Type, VOLUME_OF_FLUID
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  real,              pointer :: u_relax
  integer                    :: s, c, c1, c2
  real                       :: u_f, v_f, w_f, a12, fs
  real                       :: mass_err, dt
  real                       :: px_f, py_f, pz_f, dens_h
  character(SL)              :: solver
  real                       :: p_max, p_min, p_nor, p_nor_c
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
!     app            [m^4s/kg]
!     pp             [kg/ms^2]
!     b              [m^3/s]
!     v_flux         [m^3/s]
!
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Pressure (without solvers)')

  ! Take aliases
  grid    => flow % pnt_grid
  bulk    => flow % bulk
  v_flux  => flow % v_flux
  p       => flow % p
  pp      => flow % pp
  a       => sol % a
  b       => sol % b % val
  u_relax => flow % u_rel_corr
  dt      =  flow % dt
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, dt, ini)

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
  !-----------------------------------------!
  pp % n = 0.0

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if(c2 > 0) then

      ! Interpolate velocity
      u_f = fs * u % n(c1) + (1.0 - fs) * u % n(c2)
      v_f = fs * v % n(c1) + (1.0 - fs) * v % n(c2)
      w_f = fs * w % n(c1) + (1.0 - fs) * w % n(c2)

      ! Calculate coeficients for the system matrix
      a12 = u_relax * 0.5 * a % fc(s)                    &
                    * ( grid % vol(c1) / a % sav(c1)     &
                      + grid % vol(c2) / a % sav(c2) )

      a % val(a % pos(1,s)) = -a12
      a % val(a % pos(2,s)) = -a12
      a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) +  a12

      ! Interpolate pressure gradients as proposed by Denner
      ! (Equation 3.57 in his PhD thesis)
      dens_h = 2.0 / (1.0 / flow % density(c1) + 1.0 / flow % density(c2))
      px_f = 0.5 * dens_h * (  p % x(c1) / flow % density(c1)  &
                             + p % x(c2) / flow % density(c2) ) * grid % dx(s)
      py_f = 0.5 * dens_h * (  p % y(c1) / flow % density(c1)  &
                             + p % y(c2) / flow % density(c2) ) * grid % dy(s)
      pz_f = 0.5 * dens_h * (  p % z(c1) / flow % density(c1)  &
                             + p % z(c2) / flow % density(c2) ) * grid % dz(s)

      ! Calculate volume flux through cell face
      v_flux % n(s) = u_f * grid % sx(s)             &
                    + v_f * grid % sy(s)             &
                    + w_f * grid % sz(s)             &
                    + a12 * (p % n(c1) - p % n(c2))  &
                    + a12 * (px_f + py_f + pz_f)

      b(c1) = b(c1) - v_flux % n(s)
      b(c2) = b(c2) + v_flux % n(s)

    ! Side is on the boundary
    ! (Check: volume fluxes at the boundaries
    !  are already compute in Balance_Volume)
    else ! (c2 < 0)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        u_f = u % n(c2)
        v_f = v % n(c2)
        w_f = w % n(c2)
        v_flux % n(s) = u_f * grid % sx(s)     &
                      + v_f * grid % sy(s)     &
                      + w_f * grid % sz(s)

        b(c1) = b(c1) - v_flux % n(s)

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.   &
              Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) then

        b(c1) = b(c1) - v_flux % n(s)

        a12 = a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
        u_f = u % n(c1)
        v_f = v % n(c1)
        w_f = w % n(c1)
        v_flux % n(s) = u_f * grid % sx(s)     &
                      + v_f * grid % sy(s)     &
                      + w_f * grid % sz(s)

        b(c1) = b(c1) - v_flux % n(s)

        a12 = a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then

        b(c1) = b(c1) - v_flux % n(s)

        a12 = a % fc(s) * grid % vol(c1) / a % sav(c1)
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
    call Multiphase_Mod_Vof_Pressure_Correction(mult, sol, ini, mass_err)
  end if

  !-------------------------------------!
  !   Correct fluxes with body forces   !
  !-------------------------------------!
  call Field_Mod_Correct_Fluxes_With_Body_Forces(flow, sol)

  ! Compute mass error for SIMPLE
  if (flow % p_m_coupling == SIMPLE) then
    mass_err = 0.0
    do c = 1, grid % n_cells
      mass_err = max(mass_err, abs(b(c)))
    end do

  ! Compute mass error for PISO
  else if (flow % p_m_coupling == PISO) then
    mass_err = 0.0
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      mass_err = mass_err + abs(b(c))
    end do
    call Comm_Mod_Global_Sum_Real(mass_err)
    if (flow % i_corr == flow % n_piso_corrections) then
      flow % p % res = mass_err
      if (ini == 1) then
        flow % p % res_scal0 = mass_err
      else
        if (ini < 6) then
          flow % p % res_scal0 = max(flow % p % res_scal0, mass_err)
        end if
      end if
      flow % p % res = flow % p % res / (flow % p % res_scal0 + PICO)
    end if
  end if

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Pressure')
  if(solver .eq. 'ACM') then
    pp % tol   = PICO
    call Acm(sol,           &
             pp % n,        &
             b,             &
             pp % precond,  &
             pp % mniter,   &     ! number of V cycles
             pp % tol,      &
             pp % res,      &
             norm = p_nor)        ! last argument: number for normalisation
    stop
  else
    call Cg(sol,           &
            pp % n,        &
            b,             &
            pp % precond,  &
            pp % mniter,   &      ! max number of iterations
            pp % eniter,   &      ! executed number of iterations
            pp % tol,      &
            pp % res,      &
            norm = p_nor)         ! last argument: number for normalisation

  end if
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
  call User_Mod_End_Of_Compute_Pressure(flow, mult, dt, ini)

  call Cpu_Timer_Mod_Stop('Compute_Pressure (without solvers)')

  end subroutine
