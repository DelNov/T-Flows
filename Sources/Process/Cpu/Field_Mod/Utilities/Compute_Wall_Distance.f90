!==============================================================================!
  subroutine Compute_Wall_Distance(Flow, Sol)
!------------------------------------------------------------------------------!
!   www.cfd-online.com/Wiki/Transport_equation_based_wall_distance_calculation !
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
  integer                    :: s, c, c1, c2, n, reg
  real                       :: f_ex, f_im
  real                       :: phi_x_f, phi_y_f, phi_z_f
  real, contiguous,  pointer :: cross(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NDT = 24       ! number of false time steps
  real,    parameter :: DT  =  1.0e+6  ! false time step
!==============================================================================!

  call Profiler % Start('Compute_Wall_Distance')

  ! Take aliases
  Grid => Flow % pnt_grid
  phi  => Flow % wall_dist
  call Flow % Alias_Momentum(u, v, w)
  call Sol % Alias_Native   (A, b)

  ! If wall distance was calculated, get out of here
  if(Grid % wall_dist(1) > 0.0) then 
    call Profiler % Stop('Compute_Wall_Distance')
    return
  end if

  call Work % Connect_Real_Cell(cross)

  if(First_Proc()) then
    print '(a)',      ' # Computing wall distance ...'
    print '(a,i3,a)', ' # ... with ', NDT, ' fake time steps.'
    print '(a)',      ' # This might take a while, please wait'
  end if

  !------------------------------------------!
  !                                          !
  !                                          !
  !   Solve the equation for wall distance   !
  !                                          !
  !                                          !
  !------------------------------------------!

  ! Initialize cross diffusion term, matrix and right hand side
  cross  (:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initial values (for some cases, zero was necessary)
  do c = Cells_In_Domain()
    phi % n(c) = 0.0
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

  ! Boundary conditions for walls
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1, s)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + A % fc(s)
      end do
    end if  ! some kind of a wall
  end do    ! through boundary regions

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

    ! Source term
    do c = Cells_In_Domain()
      b(c) = b(c) + Grid % vol(c)
    end do

    ! Insert boundary conditions
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL .or.  &
         Grid % region % type(reg) .eq. WALLFL) then
        do s = Faces_In_Region(reg)
          c2 = Grid % faces_c(2, s)
          phi % n(c2) = 0.0
        end do
      else
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          phi % n(c2) = phi % n(c1)
        end do
      end if  ! is it some kind of wall?
    end do    ! through regions

    !------------------!
    !                  !
    !     Difusion     !
    !                  !
    !------------------!

    ! Gradients
    call Flow % Grad_Gauss_Variable(phi)

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

    end do  ! through faces in domain

    ! Add cross diffusion terms explicity for non-polyhedral grids
    ! (Some of them show poor convergence when solving for potential
    ! field, particularly if Grid featured concave cells near edges)
    if(.not. Grid % polyhedral) then
      do c = Cells_In_Domain()
        b(c) = b(c) + cross(c)
      end do
    end if

    !---------------------------------!
    !                                 !
    !   Solve the equations for phi   !
    !                                 !
    !---------------------------------!

    ! Call linear solver to solve the equations
    call Sol % Run(A, phi, b)

    if(First_Proc()) then
      print '(a,i4,a,e12.4)', ' # Computed wall distance in ',  phi % niter,  &
                              ' iterations with residual: ',    phi % res
    end if

    if(phi % niter .eq. 0) goto 1

    call Grid % Exchange_Cells_Real(phi % n)

    ! Re-initialize cross diffusion terms (for the next time step)
    do c = Cells_In_Domain()
      cross(c) = 0.0
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
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c2 = Grid % faces_c(2, s)
        phi % n(c2) = 0.0
      end do
    else
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        phi % n(c2) = phi % n(c1)
      end do
    end if  ! is it some kind of wall?
  end do    ! through regions

  ! Compute wall distance
  call Flow % Grad_Gauss_Variable(phi)
  do c = Cells_In_Domain_And_Buffers()
    phi % n(c) = sqrt(  phi % x(c) * phi % x(c)  &
                      + phi % y(c) * phi % y(c)  &
                      + phi % z(c) * phi % z(c)  &
                      + 2.0 * phi % n(c))        &
               - sqrt(  phi % x(c) * phi % x(c)  &
                      + phi % y(c) * phi % y(c)  &
                      + phi % z(c) * phi % z(c))
  end do

! This was used for debugging, with calculated wall distance
! call Grid % Save_Debug_Vtu(append = "wall-dist-rel-error",                   &
!                            scalar_cell = abs(phi % n - Grid % wall_dist)     &
!                                        / (Grid % wall_dist + TINY) * 100.0,  &
!                            scalar_name = "Wall Distance Error [%]")

  !---------------------------------!
  !   Save what you have computed   !
  !---------------------------------!
  Grid % wall_dist(:) = phi % n(:)

  !----------------------------------!
  !   Releasing and freeing memory   !
  !----------------------------------!
  call Work % Disconnect_Real_Cell(cross)
  call Var_Mod_Destroy_Solution(phi)

  call Profiler % Stop('Compute_Wall_Distance')

  end subroutine
