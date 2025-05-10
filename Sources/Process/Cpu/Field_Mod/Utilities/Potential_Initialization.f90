!==============================================================================!
  subroutine Potential_Initialisation(Flow, Sol, init)
!------------------------------------------------------------------------------!
!   Initializes velocity from potential (pressure-like) equation               !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  type(Solver_Type), target :: Sol
  logical, intent(in)       :: init
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: phi, u, v, w
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2, n, reg
  real                       :: f_ex, f_im
  real                       :: phi_x_f, phi_y_f, phi_z_f
  real                       :: vol_in_real, vol_in_fake, dist_min
  real, allocatable          :: store(:)
  real, contiguous,  pointer :: log_dist(:), cross(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NDT = 24       ! number of false time steps
  real,    parameter :: DT  =  1.0e+6  ! false time step
!==============================================================================!

  ! If not used, destroy the solution (release memory and solvers)
  if(.not. init) then
    call Var_Mod_Destroy_Solution(Flow % pot)
    return
  end if

  call Profiler % Start('Potential_Initialization')

  call Work % Connect_Real_Cell(log_dist, cross)

  if(First_Proc()) then
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

  ! Initialize cross diffusion term, matrix and right hand side
  cross  (:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initial values
  ! (Potential varies from 0 to 1, hence
  !  0.5 seems like a good initial guess)
  do c = Cells_In_Domain_And_Buffers()
    phi % n(c) = 0.5
  end do

  !----------------------------------------!
  !   Discretize the system of equations   !
  !----------------------------------------!

  ! Innertial term
  do c = Cells_In_Domain()
    A % val(A % dia(c)) = Grid % vol(c) / DT
  end do

  ! Diffusive term
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Fill the system matrix
    A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - A % fc(s)
    A % val(A % dia(c1))  = A % val(A % dia(c1))  + A % fc(s)
    if(Cell_In_This_Proc(c2)) then
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - A % fc(s)
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + A % fc(s)
    end if
  end do  ! through faces

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW  .or.  &
       Grid % region % type(reg) .eq. OUTFLOW .or.  &
       Grid % region % type(reg) .eq. CONVECT .or.  &
       Grid % region % type(reg) .eq. PRESSURE) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + A % fc(s)
      end do
    end if  ! region is inflow or any outlet
  end do    ! through regions

  !----------------------------------------!
  !   Begin the false time stepping loop   !
  !----------------------------------------!
  do n = 1, NDT

    ! Store the value from previous time step
    do c = Cells_In_Domain()
      phi % o(c) = phi % n(c)
    end do

    ! Innertial term
    do c = Cells_In_Domain()
      b(c) = Grid % vol(c) / DT * phi % o(c)
    end do

    ! Update boundary values
    ! (Set 1 at inflows and 0 at outflows)
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW) then
        do c2 = Cells_In_Region(reg)
          phi % n(c2) = 1.0
        end do
      else if(Grid % region % type(reg) .eq. OUTFLOW .or.  &
              Grid % region % type(reg) .eq. CONVECT .or.  &
              Grid % region % type(reg) .eq. PRESSURE) then
        do c2 = Cells_In_Region(reg)
          phi % n(c2) = 0.0
        end do
      else
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          phi % n(c2) = phi % n(c1)
        end do
      end if  ! region type
    end do    ! through regions

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
    do s = Faces_In_Domain_And_At_Buffers()

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
      cross(c1) = cross(c1) + f_ex - f_im
      cross(c2) = cross(c2) - f_ex + f_im

    end do  ! through faces

    ! Insert boundary conditions
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          b(c1) = b(c1) +  A % fc(s) * 1.0
        end do
      else if(Grid % region % type(reg) .eq. OUTFLOW .or.  &
              Grid % region % type(reg) .eq. CONVECT .or.  &
              Grid % region % type(reg) .eq. PRESSURE) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          b(c1) = b(c1) +  A % fc(s) * 0.0
        end do
      end if  ! region type
    end do    ! through regions


    ! Add cross diffusion terms explicity for non-polyhedral grids
    ! (Some of them show poor convergence when solving for potential
    ! field, particularly if Grid featured concave cells near edges)
    if(.not. Grid % polyhedral) then
      do c = Cells_In_Domain()
        b(c) = b(c) + cross(c)
      end do

      ! Fix negative sources
      do c = Cells_In_Domain()

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

    call Profiler % Start(String % First_Upper(phi % solver)  //  &
                          ' (solver for potential initialization)')

    ! Call linear solver to solve the equations
    call Sol % Run(A, phi, b)

    call Profiler % Stop(String % First_Upper(phi % solver)  //  &
                         ' (solver for potential initialization)')

    if(First_Proc()) then
      print '(a,i4,a,e12.4)', ' # Computed potential in ',   phi % niter,  &
                              ' iterations with residual: ', phi % res
    end if

    if(phi % niter .eq. 0) goto 1

    call Grid % Exchange_Cells_Real(phi % n)

    ! Recover the central coefficient in the system matrix
    if(.not. Grid % polyhedral) then
      do c = Cells_In_Domain()
        A % val(A % dia(c)) = store(c)
      end do
    end if

    ! Re-initialize cross diffusion terms (for the next time step)
    do c = Cells_In_Domain()
      cross(c) = 0.0
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
  do c = Cells_In_Domain_And_Buffers()
    log_dist(c) = Grid % wall_dist(c)
  end do
  dist_min = minval(log_dist(1:Grid % n_cells))
  call Global % Min_Real(dist_min)

  !------------------------------------------------------!
  !   Set distances from the wall to friendlier values   !
  !------------------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
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
  do c = Cells_In_Domain_And_Buffers()
    u % n(c)  = -phi % x(c) * log_dist(c)
    v % n(c)  = -phi % y(c) * log_dist(c)
    w % n(c)  = -phi % z(c) * log_dist(c)
  end do

  !-------------------------------------------------------------!
  !   Calculate real and artifical volume entering the domain   !
  !-------------------------------------------------------------!
  vol_in_real = 0.0
  vol_in_fake = 0.0
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        vol_in_real = vol_in_real + ( u % n(c2) * Grid % sx(s)    &
                                    + v % n(c2) * Grid % sy(s)    &
                                    + w % n(c2) * Grid % sz(s) )
        vol_in_fake = vol_in_fake + ( u % n(c1) * Grid % sx(s)    &
                                    + v % n(c1) * Grid % sy(s)    &
                                    + w % n(c1) * Grid % sz(s) )
      end do  ! faces at inflow
    end if    ! region is inflow
  end do      ! through regions
  call Global % Sum_Real(vol_in_real)
  call Global % Sum_Real(vol_in_fake)

  !-------------------------------------------------!
  !   Correct velocities to more realistic values   !
  !-------------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
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

  !---------------------------------!
  !   Save what you have computed   !
  !---------------------------------!
  Flow % potential(:) = phi % n(:)

  !----------------------------------!
  !   Releasing and freeing memory   !
  !----------------------------------!
  call Work % Disconnect_Real_Cell(log_dist, cross)
  call Var_Mod_Destroy_Solution(phi)

  call Profiler % Stop('Potential_Initialization')

  end subroutine
