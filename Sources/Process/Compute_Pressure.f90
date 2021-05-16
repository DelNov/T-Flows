!==============================================================================!
  subroutine Compute_Pressure(flow, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  integer, save              :: total_cells
  real                       :: total_source
  real                       :: p_max, p_min, p_nor, p_nor_c, dt, a12, fs
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
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  p      => flow % p
  pp     => flow % pp
  dt     =  flow % dt
  A      => Sol % A
  M      => Sol % M
  b      => Sol % b % val
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(flow, Vof, Sol, curr_dt, ini)

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
  A % val = 0.0

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !    (Convergence is faster with this)    !
  !-----------------------------------------!
  pp % n = 0.0

  !----------------------------------!
  !   Correct fluxes at boundaries   !
  !----------------------------------!
  call Balance_Volume(flow, Vof)

  !---------------------------------------!
  !   Compute volume fluxes at internal   !
  !    faces with Rhie and Chow method    !
  !---------------------------------------!
  call Rhie_And_Chow(flow, Vof, Sol)

  !-------------------------------------------------------------!
  !   In case of VOF, surface tension and  gravity correction   !
  !-------------------------------------------------------------!
  call Vof % Pressure_Correction(Sol)

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

        a12 = A % fc(s) * grid % vol(c1) / M % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then

        v_flux % n(s) = ( u % n(c1) * grid % sx(s)     &
                        + v % n(c1) * grid % sy(s)     &
                        + w % n(c1) * grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

        a12 = A % fc(s) * grid % vol(c1) / M % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12

      else  ! it is SYMMETRY
        v_flux % n(s) = 0.0
      end if
    end if

  end do

  !-------------------------------------------------------------------------!
  !   In case of mass transfer, add addtional source to pressure equation   !
  !-------------------------------------------------------------------------!
  call Vof % Mass_Transfer_Pressure_Source(b)

  !----------------------------------------!
  !   Balance the source over processors   !
  !----------------------------------------!
  ! Since the introduction of face-based body forces, a very slight
  ! imbalance in pressure source can occur globally (over all processors)
  ! and we are talking about the values of the order 1e-18 here.  This
  ! happens because of the round-off errors associated with real arithmetics
  ! performed in CPUs, meaning that: A + B = B + A + epsilon.  If forces
  ! are computed over faces, they are slightly different in different
  ! processors due to these round-off errors, although values in cells
  ! surrounding them are the same after exchanging the buffers.  In order
  ! to fix it, the balancing procedure which follows is introduced.
  ! However, the fix should not be applied if domain has pressure outlet.
  if( .not. flow % has_pressure_outlet) then
    if(total_cells .eq. 0) then  ! wasn't set yet
      total_cells = grid % n_cells - grid % comm % n_buff_cells
      call Comm_Mod_Global_Sum_Int(total_cells)
    end if
    total_source = sum(b(1:(grid % n_cells - grid % comm % n_buff_cells)))
    call Comm_Mod_Global_Sum_Real(total_source)
    b(:) = b(:) - total_source/real(total_cells)
  end if

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer % Start('Linear_Solver_For_Pressure')
  call Sol % Cg(A,             &
                pp % n,        &
                b,             &
                pp % precond,  &
                pp % mniter,   &      ! max number of iterations
                pp % eniter,   &      ! executed number of iterations
                pp % tol,      &
                pp % res,      &
                norm = p_nor)         ! number for normalisation

  call Cpu_Timer % Stop('Linear_Solver_For_Pressure')

  if (flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
  else
    if (flow % i_corr == flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
    end if
  end if

  call Field_Mod_Grad_Pressure_Correction(flow, pp)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = -grid % n_bnd_cells, grid % n_cells
    p % n(c) =  p % n(c) + pp % urf * pp % n(c)
  end do

  call Field_Mod_Grad_Pressure(flow, p)

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:grid % n_cells))
  p_min  = minval(p % n(1:grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(flow, Vof, Sol, curr_dt, ini)

  call Cpu_Timer % Stop('Compute_Pressure (without solvers)')

  end subroutine
