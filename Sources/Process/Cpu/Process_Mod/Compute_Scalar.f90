!==============================================================================!
  subroutine Compute_Scalar(Process, Flow, Turb, Vof, Sol, sc)
!------------------------------------------------------------------------------!
!>  The subroutine Compute_Scalar in T-Flows is dedicated to solving transport
!>  equations for scalar quantities such as concentrations or passive scalars.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and setup: Begins by setting up necessary local variables !
!     and pointers to important data structures like the grid and scalar field.!
!     This includes establishing aliases for efficient data access and         !
!     starting a profiler for performance monitoring.                          !
!   * Equation discretization: It discretizes the scalar transport equations,  !
!     computing terms related to advection and diffusion. This includes        !
!     handling gradients of scalar quantities and the computation of diffusive !
!     fluxes.                                                                  !
!   * Spatial discretization and matrix assembly: Iterates through all grid    !
!     faces to calculate coefficients of the discretized scalar equation.      !
!     Assembles the system matrix for the scalar equations, considering        !
!     diffusive fluxes and cross diffusion.                                    !
!   * Under-Relaxation and solving equations: Applies under-relaxation to      !
!     the equations and then uses the linear solver to solve the discretized   !
!     equations for the scalar quantities. The subroutine handles both         !
!     explicit and implicit parts of diffusive fluxes.                         !
!   * Post-processing and boundary updates: After solving the equations,       !
!     updates the boundary values for scalar quantities and refreshes buffers. !
!     This step ensures that boundary conditions are correctly applied to the  !
!     newly computed scalar fields. Also includes gradient computation         !
!     post-solution.                                                           !
!   * User-defined functions: Throughout its execution, Compute_Scalar invokes !
!     various user-defined functions for additional customization. These       !
!     functions enable users to integrate specific behaviors or calculations   !
!     at different stages of the scalar computation process.                   !
!   * Performance monitoring: The subroutine consistently monitors its         !
!     performance, contributing to the analysis and optimization of the        !
!     simulation. This monitoring helps in identifying bottlenecks and         !
!     optimizing the simulation process for better efficiency.                 !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Turb_Type),     target :: Turb     !! turbulence object
  type(Vof_Type),      target :: Vof      !! VOF object
  type(Solver_Type),   target :: Sol      !! solver object
  integer, intent(in)         :: sc       !! scalar rank
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: phi
  integer                    :: c, s, c1, c2, row, col
  real                       :: a12, a21
  real                       :: rs, dt
  real                       :: dif_eff, f_ex, f_im
  real                       :: phi_stress, q_exp
  real                       :: phix_f, phiy_f, phiz_f
  real, contiguous,  pointer :: q_turb(:), cross(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                /
!    |     d phi      |                |
!    | rho ----- dV   | rho u phi dS = | gamma DIV phi dS
!    |      dt        |                |
!   /                /                /
!
!==============================================================================!

  call Profiler % Start('Compute_Scalar (without solvers)')

  call Work % Connect_Real_Cell(q_turb, cross)

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  phi    => Flow % scalar(sc)
  dt     =  Flow % dt
  call Turb % Alias_Stresses(uu, vv, ww, uv, uw, vw)
  call Sol % Alias_Native   (A, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Scalar(Flow, Turb, Vof, Sol, sc)

  ! Initialize cross diffusion sources, matrix and right hand side
  cross(:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initialize turbulent scalar fluxes
  ! and fluxes coming from interfaces
  q_turb(:) = 0.0

  !--------------------------!
  !   Initialize variables   !
  !--------------------------!

  ! Old values (o and oo)
  if(Iter % Current() .lt. 2) then
    do c = Cells_In_Domain_And_Buffers()
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Flow % Grad_Variable(phi)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, Flow % density, v_flux % n, b)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!

  ! It used to read sc_t from here which is an overkill, so check
  Assert(sc_t > 0.0)

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    call Turb % Face_Diff_And_Stress(dif_eff, phi_stress, s, sc)

    ! Gradients on the cell face
    phix_f = Grid % fw(s)*phi % x(c1) + (1.0-Grid % fw(s))*phi % x(c2)
    phiy_f = Grid % fw(s)*phi % y(c1) + (1.0-Grid % fw(s))*phi % y(c2)
    phiz_f = Grid % fw(s)*phi % z(c1) + (1.0-Grid % fw(s))*phi % z(c2)

    ! Total (exact) diffusive flux
    f_ex = dif_eff  * (  phix_f * Grid % sx(s)  &
                       + phiy_f * Grid % sy(s)  &
                       + phiz_f * Grid % sz(s))

    ! Implicit diffusive flux
    f_im = dif_eff  * A % fc(s)          &
         * (  phix_f * Grid % dx(s)      &
            + phiy_f * Grid % dy(s)      &
            + phiz_f * Grid % dz(s) )

    ! Calculate the coefficients for the sysytem matrix
    a12 = dif_eff * A % fc(s)
    a21 = dif_eff * A % fc(s)

    ! Blend system matrix if desired to do so
    if(phi % blend_matrix) then
      a12 = a12  - min(v_flux % n(s), 0.0) * Flow % density(c1)
      a21 = a21  + max(v_flux % n(s), 0.0) * Flow % density(c2)
    end if

    ! Cross diffusion part
    cross(c1) = cross(c1) + f_ex - f_im
    if(c2 .gt. 0) then
      cross(c2) = cross(c2) - f_ex + f_im
    end if

    ! Put the influence of turbulent scalar fluxes explicitly in the system
    q_turb(c1) = q_turb(c1) + phi_stress
    if(c2 > 0) then
      q_turb(c2) = q_turb(c2) - phi_stress
    end if

    ! Fill the system matrix
    if(c2 > 0) then
      A % val(A % dia(c1))  = A % val(A % dia(c1)) + a12
      A % val(A % dia(c2))  = A % val(A % dia(c2)) + a21
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
    else if(c2 < 0) then

      ! Conditions which are always of Dirichlet type
      if( (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. CONVECT) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * phi % n(c2)

      ! Ambient when it is inflow (see the v_flux check)
      else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. AMBIENT  &
              .and. v_flux % n(s) .lt. 0.0) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1)  = b(c1) + a12 * phi % n(c2)  ! phi % n(c2) is ambient value here

      else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
        b(c1) = b(c1) + Grid % s(s) * phi % q(c2)

      end if  ! boundary condition type
    end if    ! c2 .lt. 0
  end do      ! through faces

  !----------------------------------------------!
  !   Explicitly treated diffusion scalar fluxes !
  !   and cross diffusion                        !
  !----------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()

    ! Total explicit heat flux
    q_exp = cross(c) + q_turb(c)

    if(q_exp >= 0) then
      b(c)  = b(c) + q_exp
    else
      A % val(A % dia(c)) = A % val(A % dia(c)) - q_exp / (phi % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Numerics_Mod_Inertial_Term(phi, Flow % density, A, b, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!

  call User_Mod_Source(Flow, phi, A, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, A, b)

  call Profiler % Start(String % First_Upper(phi % solver)  //  &
                        ' (solver for scalars)')

  ! Call linear solver to solve them
  call Sol % Run(A, phi, b)

  call Profiler % Stop(String % First_Upper(phi % solver)  //  &
                       ' (solver for scalars)')

  rs = sc                      ! reterive the rank of scalar
  row = ceiling(rs/6)          ! will be 1 (scal. 1-6), 2 (scal. 6-12), etc.
  col = nint(rs) - (row-1)*6   ! will be in range 1 - 6

  call Info % Iter_Fill_Scalar_At(row, col, phi % name, phi % res, phi % niter)

  call Flow % Grad_Variable(phi)

  ! User function
  call User_Mod_End_Of_Compute_Scalar(Flow, Turb, Vof, Sol, sc)

  call Work % Disconnect_Real_Cell(q_turb, cross)

  call Profiler % Stop('Compute_Scalar (without solvers)')

  end subroutine
