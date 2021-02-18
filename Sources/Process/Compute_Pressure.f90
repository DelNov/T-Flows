!==============================================================================!
  subroutine Compute_Pressure(flow, mult, sol, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: b_save => r_cell_01
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
  real                       :: u_f, v_f, w_f, a12, fs, dt
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
!     app              [m^4s/kg]
!     pp               [kg/ms^2]
!     p%x, p%y, p%z    [kg/m^2s^2]
!     px_f, py_f, pz_f [kg/m^2s^2]
!     b                [m^3/s]
!     v_flux           [m^3/s]
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
  call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, ini)

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

      ! If there is a jump in velocities, call specialized gradient calculation
      if(mult % mass_transfer) then
        u_f = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, u, s)
        v_f = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, v, s)
        w_f = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, w, s)

      ! No jumps, call usual routines
      else
        u_f = Field_Mod_Interpolate_Var_To_Face(flow, u, s)
        v_f = Field_Mod_Interpolate_Var_To_Face(flow, v, s)
        w_f = Field_Mod_Interpolate_Var_To_Face(flow, w, s)
      end if

      ! Calculate coeficients for the system matrix
      ! a12 [m*m^3*s/kg = m^4s/kg]
      if(.not. mult % mass_transfer) then
        a12 = u_relax * 0.5 * a % fc(s)                    &
                      * ( grid % vol(c1) / a % sav(c1)     &
                        + grid % vol(c2) / a % sav(c2) )
      else
        if(mult % cell_at_elem(c1) .eq. 0 .and.  &
           mult % cell_at_elem(c2) .eq. 0 .or.   &
           mult % cell_at_elem(c1) .ne. 0 .and.  &
           mult % cell_at_elem(c2) .ne. 0) then
          a12 = u_relax * 0.5 * a % fc(s)                    &
                        * ( grid % vol(c1) / a % sav(c1)     &
                          + grid % vol(c2) / a % sav(c2) )
        else
          if(mult % cell_at_elem(c1) .eq. 0 .and.  &
             mult % cell_at_elem(c2) .ne. 0) then
            a12 = u_relax * a % fc(s) * grid % vol(c1) / a % sav(c1)
          end if
          if(mult % cell_at_elem(c1) .ne. 0 .and.  &
             mult % cell_at_elem(c2) .eq. 0) then
            a12 = u_relax * a % fc(s) * grid % vol(c2) / a % sav(c2)
          end if
        end if
      end if

      a % val(a % pos(1,s)) = -a12
      a % val(a % pos(2,s)) = -a12
      a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) +  a12

      ! Interpolate pressure gradients as proposed by Denner
      ! (Equation 3.57 in his PhD thesis)
      ! dens_h           [kg/m^3]
      ! px_f, py_f, pz_f [kg/m^2s^2]
      dens_h = 2.0 / (1.0 / flow % density(c1) + 1.0 / flow % density(c2))
      px_f = 0.5 * dens_h * (  p % x(c1) / flow % density(c1)  &
                             + p % x(c2) / flow % density(c2) )
      py_f = 0.5 * dens_h * (  p % y(c1) / flow % density(c1)  &
                             + p % y(c2) / flow % density(c2) )
      pz_f = 0.5 * dens_h * (  p % z(c1) / flow % density(c1)  &
                             + p % z(c2) / flow % density(c2) )

      ! Calculate current volume flux through cell face with pressure
      ! defined at a cell face and assuming that pressure correction
      ! (pp) part is treated implicitly
      v_flux % n(s) = u_f * grid % sx(s)             &
                    + v_f * grid % sy(s)             &
                    + w_f * grid % sz(s)             &
                    + a12 * (p % n(c1) - p % n(c2))  &
                    + a12 * (  px_f * grid % dx(s)   &
                             + py_f * grid % dy(s)   &
                             + pz_f * grid % dz(s))

      ! Any of the cells is at interface, use only non-interface value
      ! (This is the old way, and seems to be working better after all)
      if(mult % mass_transfer) then
        if(mult % cell_at_elem(c1) .eq. 0 .and.  &
           mult % cell_at_elem(c2) .ne. 0) then
          v_flux % n(s) = u_f * grid % sx(s)     &
                        + v_f * grid % sy(s)     &
                        + w_f * grid % sz(s)
        end if
        if(mult % cell_at_elem(c1) .ne. 0 .and.  &
           mult % cell_at_elem(c2) .eq. 0) then
          v_flux % n(s) = u_f * grid % sx(s)     &
                        + v_f * grid % sy(s)     &
                        + w_f * grid % sz(s)
        end if
      end if

      b(c1) = b(c1) - v_flux % n(s)
      b(c2) = b(c2) + v_flux % n(s)

    ! Face is on the boundary
    ! (Check: volume fluxes at the boundaries was
    !  corrected in Balance_Volume called above)
    else ! (c2 < 0)

      b(c1) = b(c1) - v_flux % n(s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.   &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) then

        a12 = a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then

        a12 = a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

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
  call Field_Mod_Correct_Fluxes_With_Body_Forces(flow, sol)

  ! Compute volume error
  flow % vol_res = 0.0
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    flow % vol_res = flow % vol_res + abs(b(c))
  end do
  call Comm_Mod_Global_Sum_Real(flow % vol_res)

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Pressure')
  call Solver_Mod_Cg(sol,           &
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
  call User_Mod_End_Of_Compute_Pressure(flow, mult, ini)

  call Cpu_Timer_Mod_Stop('Compute_Pressure (without solvers)')

  end subroutine
