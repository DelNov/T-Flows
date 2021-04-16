!==============================================================================!
  subroutine Field_Mod_Potential_Initialization(flow, sol)
!------------------------------------------------------------------------------!
!   Discretizes and solves eliptic relaxation equations for f22.               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: log_dist => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: phi, u, v, w
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2, n
  real                       :: f_ex, f_im
  real                       :: phi_x_f, phi_y_f, phi_z_f
  real                       :: vol_in_real, vol_in_fake, dist_min
  real, allocatable          :: store(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NDT = 24       ! number of false time steps
  real,    parameter :: DT  =  1.0e+6  ! false time step
!==============================================================================!

  if(this_proc < 2) then
    print '(a)',      ' # Computing potential to initialize velocity field ...'
    print '(a,i3,a)', ' # ... with ', NDT, ' fake time steps.'
    print '(a)',      ' # This might take a while, please wait'
  end if

  ! Take aliases
  grid => flow % pnt_grid
  phi  => flow % pot

  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Solver_Mod_Alias_System (sol,  a, b)
  allocate(store(grid % n_cells)); store(:) = 0.0

  !--------------------------------------!
  !                                      !
  !                                      !
  !   Solve the equation for potential   !
  !                                      !
  !                                      !
  !--------------------------------------!

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! Initial values
  ! (Potential varies from 0 to 1, hence
  !  0.5 seems like a good initial guess)
  do c = 1, grid % n_cells
    phi % n(c) = 0.5
  end do

  !----------------------------------------!
  !   Discretize the system of equations   !
  !----------------------------------------!

  ! Innertial term
  do c = 1, grid % n_cells
    a % val(a % dia(c)) = grid % vol(c) / DT
  end do

  ! Diffusive term
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Fill the system matrix
    if(c2  > 0) then
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a % fc(s)
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a % fc(s)
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a % fc(s)
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a % fc(s)
    else if(c2  < 0) then
      ! Inflow and outflow
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)   .or.  &
          (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) .or.  &
          (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
          (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE)) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a % fc(s)
      end if
    end if
  end do  ! through faces

  !----------------------------------------!
  !   Begin the false time stepping loop   !
  !----------------------------------------!
  do n = 1, NDT

    ! Store the value from previous time step
    do c = 1, grid % n_cells
      phi % o(c) = phi % n(c)
    end do

    ! Innertial term
    do c = 1, grid % n_cells
      b(c) = grid % vol(c) / DT * phi % o(c)
    end do

    ! Update boundary values
    ! (Set 1 at inflows and 0 at outflows)
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)

      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)) then
          phi % n(c2) = 1.0
        else if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) .or.  &
                 (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
                 (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE)) then
          phi % n(c2) = 0.0
        else
          phi % n(c2) = phi % n(c1)
        end if
      end if
    end do

    !------------------!
    !                  !
    !     Difusion     !
    !                  !
    !------------------!

    ! Gradients
    call Field_Mod_Grad_Variable(flow, phi)

    !----------------------------!
    !   Spatial discretization   !
    !----------------------------!
    do s = 1, grid % n_faces

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      phi_x_f = grid % fw(s) * phi % x(c1) + (1.0-grid % fw(s)) * phi % x(c2)
      phi_y_f = grid % fw(s) * phi % y(c1) + (1.0-grid % fw(s)) * phi % y(c2)
      phi_z_f = grid % fw(s) * phi % z(c1) + (1.0-grid % fw(s)) * phi % z(c2)

      ! Total (exact) diffusive flux
      f_ex = (  phi_x_f * grid % sx(s)   &
              + phi_y_f * grid % sy(s)   &
              + phi_z_f * grid % sz(s) )

      ! Implicit diffusive flux
      f_im=(   phi_x_f * grid % dx(s)        &
             + phi_y_f * grid % dy(s)        &
             + phi_z_f * grid % dz(s)) * a % fc(s)

      ! Cross diffusion part
      phi % c(c1) = phi % c(c1) + f_ex - f_im
      if(c2  > 0) then
        phi % c(c2) = phi % c(c2) - f_ex + f_im
      end if

      ! Fill the system matrix
      if(c2 < 0) then

        ! Inflow
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)) then
          b(c1) = b(c1) +  a % fc(s) * 1.0
        end if

        ! Outflow
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) .or.  &
            (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
            (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE)) then
          b(c1) = b(c1) + a % fc(s) * 0.0
        end if

        ! For wall and wall flux, solid walls in any case => do nothing!
      end if

    end do  ! through faces

    ! Add cross diffusion terms explicity for non-polyhedral grids
    ! (Some of them show poor convergence when solving for potential
    ! field, particularly if grid featured concave cells near edges)
    if(.not. grid % polyhedral) then
      do c = 1, grid % n_cells
        b(c) = b(c) + phi % c(c)
      end do

      ! Fix negative sources
      do c = 1, grid % n_cells

        ! Store the central term
        store(c) = a % val(a % dia(c))

        ! If source is negative, fix it!
        if(b(c) < 0.0) then
          a % val(a % dia(c)) = a % val(a % dia(c)) - b(c) / phi % o(c)
          b(c) = 0.0
        end if
      end do
    end if

    !---------------------------------!
    !                                 !
    !   Solve the equations for phi   !
    !                                 !
    !---------------------------------!

    ! Set number of iterations "by hand"
    phi % mniter  = 99
    phi % tol     =  1.0e-5
    phi % precond = 'INCOMPLETE_CHOLESKY'

    ! Call linear solver to solve the equations
    call Solver_Mod_Bicg(sol,            &
                         a,              &
                         phi % n,        &
                         b,              &
                         phi % precond,  &
                         phi % mniter,   &
                         phi % eniter,   &
                         phi % tol,      &
                         phi % res)
    if(this_proc < 2) then
      print '(a,i4,a,e12.4)', ' # Computed potential in ',   phi % eniter,  &
                              ' iterations with residual: ', phi % res
    end if

    if(phi % eniter .eq. 0) goto 1

    call Grid_Mod_Exchange_Cells_Real(grid, phi % n)

    ! Recover the central coefficient in the system matrix
    if(.not. grid % polyhedral) then
      do c = 1, grid % n_cells
        a % val(a % dia(c)) = store(c)
      end do
    end if

    ! Re-initialize cross diffusion terms (for the next time step)
    do c = 1, grid % n_cells
      phi % c(c) = 0.0
    end do

  end do

  !-------------------------------!
  !                               !
  !                               !
  !   Compute modified distance   !
  !                               !
  !                               !
  !-------------------------------!
1 continue

  !---------------------------------------------!
  !   Find the minimum distance from the wall   !
  !---------------------------------------------!
  do c = 1, grid % n_cells
    log_dist(c) = grid % wall_dist(c)
  end do
  dist_min = minval(log_dist(1:grid % n_cells))
  call Comm_Mod_Global_Min_Real(dist_min)

  !------------------------------------------------------!
  !   Set distances from the wall to friendlier values   !
  !------------------------------------------------------!
  do c = 1, grid % n_cells
    log_dist(c) = log( log_dist(c) / (0.99 * dist_min) )
  end do

  !--------------------------------!
  !                                !
  !                                !
  !   Set initial velocity field   !
  !                                !
  !                                !
  !--------------------------------!

  call Field_Mod_Grad_Variable(flow, phi)

  !---------------------------------------!
  !   Set first estimate for velocities   !
  !---------------------------------------!
  do c = 1, grid % n_cells
    u % n(c)  = -phi % x(c) * log_dist(c)
    v % n(c)  = -phi % y(c) * log_dist(c)
    w % n(c)  = -phi % z(c) * log_dist(c)
  end do

  !-------------------------------------------------------------!
  !   Calculate real and artifical volume entering the domain   !
  !-------------------------------------------------------------!
  vol_in_real = 0.0
  vol_in_fake = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. INFLOW) then
        vol_in_real = vol_in_real + ( u % n(c2)*grid % sx(s)    &
                                    + v % n(c2)*grid % sy(s)    &
                                    + w % n(c2)*grid % sz(s) )
        vol_in_fake = vol_in_fake + ( u % n(c1)*grid % sx(s)    &
                                    + v % n(c1)*grid % sy(s)    &
                                    + w % n(c1)*grid % sz(s) )
      end if
    end if
  end do
  call Comm_Mod_Global_Sum_Real(vol_in_real)
  call Comm_Mod_Global_Sum_Real(vol_in_fake)

  !-------------------------------------------------!
  !   Correct velocities to more realistic values   !
  !-------------------------------------------------!
  do c = 1, grid % n_cells
    u % n(c)  = u % n(c) * vol_in_real / vol_in_fake
    v % n(c)  = v % n(c) * vol_in_real / vol_in_fake
    w % n(c)  = w % n(c) * vol_in_real / vol_in_fake
    u % o(c)  = u % n(c)
    v % o(c)  = v % n(c)
    w % o(c)  = w % n(c)
    u % oo(c) = u % n(c)
    v % oo(c) = v % n(c)
    w % oo(c) = w % n(c)
  end do

  end subroutine
