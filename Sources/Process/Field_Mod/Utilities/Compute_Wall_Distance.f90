!==============================================================================!
  subroutine Compute_Wall_Distance(Flow, Sol)
!------------------------------------------------------------------------------!
!   www.cfd-online.com/Wiki/Transport_equation_based_wall_distance_calculation !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
! use Work_Mod, only: log_dist => r_cell_01
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                     (allocates Work_Mod)          !
!     |                                                                        !
!     +----> Field_Mod_Potential_Initialization  (safe to use r_cell_01)       !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  type(Solver_Type), target :: Sol
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: phi, u, v, w
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2, n
  real                       :: f_ex, f_im
  real                       :: phi_x_f, phi_y_f, phi_z_f
  real                       :: vol_in_real, vol_in_fake, dist_min
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
  Grid => Flow % pnt_grid
  phi  => Flow % wall_dist

  call Flow % Alias_Momentum(u, v, w)
  call Sol % Alias_Native   (A, b)

  !------------------------------------------!
  !                                          !
  !                                          !
  !   Solve the equation for wall distance   !
  !                                          !
  !                                          !
  !------------------------------------------!

  ! Initialize matrix and right hand side
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initial values
  ! (Potential varies from 0 to 1, hence
  !  0.5 seems like a good initial guess)
  do c = 1, Grid % n_cells
    phi % n(c) = 0.5
  end do

  !----------------------------------------!
  !   Discretize the system of equations   !
  !----------------------------------------!

  ! Innertial term
  do c = 1, Grid % n_cells
    A % val(A % dia(c)) = Grid % vol(c) / DT
  end do

  ! Diffusive term
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Fill the system matrix
    if(c2  > 0) then
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A % fc(s)
      A % val(A % dia(c1))  = A % val(A % dia(c1))  + A % fc(s)
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A % fc(s)
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + A % fc(s)
    else if(c2  < 0) then
      ! Walls
      if( (Grid % Bnd_Cond_Type(c2) .eq. WALL)   .or.  &
          (Grid % Bnd_Cond_Type(c2) .eq. WALLFL)) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + A % fc(s)
      end if
    end if
  end do  ! through faces

  !----------------------------------------!
  !   Begin the false time stepping loop   !
  !----------------------------------------!
  do n = 1, NDT

    ! Store the value from previous time step
    do c = 1, Grid % n_cells
      phi % o(c) = phi % n(c)
    end do

    ! Innertial term
    do c = 1, Grid % n_cells
      b(c) = Grid % vol(c) / DT * phi % o(c)
    end do

    ! Source term
    do c = 1, Grid % n_cells
      b(c) = b(c) + Grid % vol(c)
    end do

    ! Update boundary values
    ! (Set 1 at inflows and 0 at outflows)
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      if(c2 < 0) then
        ! On all walls set to zero
        if( (Grid % Bnd_Cond_Type(c2) .eq. WALL)   .or.  &
            (Grid % Bnd_Cond_Type(c2) .eq. WALLFL)) then
          phi % n(c2) = 0.0

        ! Copy from inside everywhere else
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
    call Flow % Grad_Variable(phi)

    !----------------------------!
    !   Spatial discretization   !
    !----------------------------!
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      phi_x_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
      phi_y_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
      phi_z_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

      ! Total (exact) diffusive flux
      f_ex = (  phi_x_f * Grid % sx(s)   &
              + phi_y_f * Grid % sy(s)   &
              + phi_z_f * Grid % sz(s) )

      ! Implicit diffusive flux
      f_im=(   phi_x_f * Grid % dx(s)        &
             + phi_y_f * Grid % dy(s)        &
             + phi_z_f * Grid % dz(s)) * A % fc(s)

      ! Cross diffusion part
      phi % c(c1) = phi % c(c1) + f_ex - f_im
      if(c2  > 0) then
        phi % c(c2) = phi % c(c2) - f_ex + f_im
      end if

      ! Fill the system matrix
      if(c2 < 0) then
        ! For wall and wall flux, solid walls in any case => do nothing!
      end if

    end do  ! through faces

    ! Add cross diffusion terms explicity for non-polyhedral grids
    ! (Some of them show poor convergence when solving for potential
    ! field, particularly if Grid featured concave cells near edges)
    if(.not. Grid % polyhedral) then
      do c = 1, Grid % n_cells
        b(c) = b(c) + phi % c(c)
      end do
    end if

    !---------------------------------!
    !                                 !
    !   Solve the equations for phi   !
    !                                 !
    !---------------------------------!

    ! Call linear solver to solve the equations
    call Sol % Run(phi % solver,     &
                   phi % prec,       &
                   phi % prec_opts,  &
                   A,                &
                   phi % n,          &
                   b,                &
                   phi % mniter,     &
                   phi % eniter,     &
                   phi % tol,        &
                   phi % res)

    if(this_proc < 2) then
      print '(a,i4,a,e12.4)', ' # Computed wall distance in ',  phi % eniter,  &
                              ' iterations with residual: ',    phi % res
    end if

    if(phi % eniter .eq. 0) goto 1

    call Grid % Exchange_Cells_Real(phi % n)

    ! Re-initialize cross diffusion terms (for the next time step)
    do c = 1, Grid % n_cells
      phi % c(c) = 0.0
    end do

  end do

  !---------------------------!
  !                           !
  !                           !
  !   Compute wall distance   !
  !                           !
  !                           !
  !---------------------------!
1 continue

  ! Set it to zero on boundaries (probably not needed)
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 < 0) then
      ! On all walls set to zero
      if( (Grid % Bnd_Cond_Type(c2) .eq. WALL)   .or.  &
          (Grid % Bnd_Cond_Type(c2) .eq. WALLFL)) then
        phi % n(c2) = 0.0

      ! Copy from inside everywhere else
      else
        phi % n(c2) = phi % n(c1)
      end if
    end if
  end do

  ! Compute wall distance
  call Flow % Grad_Variable(phi)
  do c = 1, Grid % n_cells
    phi % n(c) = sqrt(  phi % x(c) * phi % x(c)  &
                      + phi % y(c) * phi % y(c)  &
                      + phi % z(c) * phi % z(c)  &
                      + 2.0 * phi % n(c))        &
               - sqrt(  phi % x(c) * phi % x(c)  &
                      + phi % y(c) * phi % y(c)  &
                      + phi % z(c) * phi % z(c))
  end do
  call Save_Debug_Vtu(Grid, "wall_dist_equation",          &
                             scalar_cell = phi % n,        &
                             scalar_name = "wall_dist_equation")
  call Save_Debug_Vtu(Grid, "wall_dist_error",                              &
                             scalar_cell = abs(phi % n - Grid % wall_dist)  &
                                         / Grid % wall_dist * 100.0,        &
                             scalar_name = "wall_dist_error [%]")

  end subroutine
