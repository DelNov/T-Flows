!==============================================================================!
  subroutine Potential_Initialization(Flow, Sol)
!------------------------------------------------------------------------------!
!   Initializes velocity from potential (pressure-like) equation               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: log_dist => r_cell_01
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
  Grid => Flow % pnt_grid
  phi  => Flow % pot

  call Flow % Alias_Momentum(u, v, w)
  call Sol % Alias_Native(A, b)
  allocate(store(Grid % n_cells)); store(:) = 0.0

  !--------------------------------------!
  !                                      !
  !                                      !
  !   Solve the equation for potential   !
  !                                      !
  !                                      !
  !--------------------------------------!

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
      ! Inflow and outflow
      if( (Grid % Bnd_Cond_Type(c2) .eq. INFLOW)   .or.  &
          (Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) .or.  &
          (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) .or.  &
          (Grid % Bnd_Cond_Type(c2) .eq. PRESSURE)) then
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

    ! Update boundary values
    ! (Set 1 at inflows and 0 at outflows)
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. INFLOW)) then
          phi % n(c2) = 1.0
        else if( (Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) .or.  &
                 (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) .or.  &
                 (Grid % Bnd_Cond_Type(c2) .eq. PRESSURE)) then
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

        ! Inflow
        if( (Grid % Bnd_Cond_Type(c2) .eq. INFLOW)) then
          b(c1) = b(c1) +  A % fc(s) * 1.0
        end if

        ! Outflow
        if( (Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) .or.  &
            (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) .or.  &
            (Grid % Bnd_Cond_Type(c2) .eq. PRESSURE)) then
          b(c1) = b(c1) + A % fc(s) * 0.0
        end if

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

      ! Fix negative sources
      do c = 1, Grid % n_cells

        ! Store the central term
        store(c) = A % val(A % dia(c))

        ! If source is negative, fix it!
        if(b(c) < 0.0) then
          A % val(A % dia(c)) = A % val(A % dia(c)) - b(c) / phi % o(c)
          b(c) = 0.0
        end if
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
      print '(a,i4,a,e12.4)', ' # Computed potential in ',   phi % eniter,  &
                              ' iterations with residual: ', phi % res
    end if

    if(phi % eniter .eq. 0) goto 1

    call Grid % Exchange_Cells_Real(phi % n)

    ! Recover the central coefficient in the system matrix
    if(.not. Grid % polyhedral) then
      do c = 1, Grid % n_cells
        A % val(A % dia(c)) = store(c)
      end do
    end if

    ! Re-initialize cross diffusion terms (for the next time step)
    do c = 1, Grid % n_cells
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
  do c = 1, Grid % n_cells
    log_dist(c) = Grid % wall_dist(c)
  end do
  dist_min = minval(log_dist(1:Grid % n_cells))
  call Comm_Mod_Global_Min_Real(dist_min)

  !------------------------------------------------------!
  !   Set distances from the wall to friendlier values   !
  !------------------------------------------------------!
  do c = 1, Grid % n_cells
    log_dist(c) = log( log_dist(c) / (0.99 * dist_min) )
  end do

  !--------------------------------!
  !                                !
  !                                !
  !   Set initial velocity field   !
  !                                !
  !                                !
  !--------------------------------!

  call Flow % Grad_Variable(phi)

  !---------------------------------------!
  !   Set first estimate for velocities   !
  !---------------------------------------!
  do c = 1, Grid % n_cells
    u % n(c)  = -phi % x(c) * log_dist(c)
    v % n(c)  = -phi % y(c) * log_dist(c)
    w % n(c)  = -phi % z(c) * log_dist(c)
  end do

  !-------------------------------------------------------------!
  !   Calculate real and artifical volume entering the domain   !
  !-------------------------------------------------------------!
  vol_in_real = 0.0
  vol_in_fake = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2  < 0) then
      if(Grid % Bnd_Cond_Type( c2) .eq. INFLOW) then
        vol_in_real = vol_in_real + ( u % n(c2)*Grid % sx(s)    &
                                    + v % n(c2)*Grid % sy(s)    &
                                    + w % n(c2)*Grid % sz(s) )
        vol_in_fake = vol_in_fake + ( u % n(c1)*Grid % sx(s)    &
                                    + v % n(c1)*Grid % sy(s)    &
                                    + w % n(c1)*Grid % sz(s) )
      end if
    end if
  end do
  call Comm_Mod_Global_Sum_Real(vol_in_real)
  call Comm_Mod_Global_Sum_Real(vol_in_fake)

  !-------------------------------------------------!
  !   Correct velocities to more realistic values   !
  !-------------------------------------------------!
  do c = 1, Grid % n_cells
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
