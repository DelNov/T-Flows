!==============================================================================!
  subroutine Compute_Energy(flow, sol, dt, ini)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Control_Mod
  use Field_Mod,    only: Field_Type, conductivity, capacity, density
  use Rans_Mod
  use Var_Mod,      only: Var_Type
  use Grid_Mod,     only: Grid_Type
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod, only: CENTRAL, LINEAR, PARABOLIC
  use Solver_Mod,   only: Solver_Type, Bicg, Cg, Cgs
  use Matrix_Mod,   only: Matrix_Type
  use Work_Mod,     only: t_min => r_cell_04,  &
                          t_max => r_cell_05
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  integer                   :: ini
  real                      :: dt
!----------------------------------[Calling]-----------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------! 
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer           :: n, c, s, c1, c2, niter
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: con_eff1, f_ex1, f_im1, tx_f1, ty_f1, tz_f1
  real              :: con_eff2, f_ex2, f_im2, tx_f2, ty_f2, tz_f2
  real              :: ts, pr_t1, pr_t2
  real              :: ut_s, vt_s, wt_s, t_stress, con_t
  character(len=80) :: precond       ! preconditioner
  integer           :: adv_scheme    ! space-discretiztion of advection scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  integer           :: td_scheme     ! time-disretization for inerita  
  real              :: urf           ! under-relaxation factor
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                 /
!    |        dT      |                 |
!    | rho Cp -- dV   | rho u Cp T dS = | lambda  DIV T dS
!    |        dt      |                 |
!   /                /                 /
!
!
!  Dimension of the system under consideration
!
!     [A]{T} = {b}   [J/s = W]  
!
!  Dimensions of certain variables:
!
!     Cp     [J/kg K] (heat capacity)
!     lambda [W/m K] (heat conductivity)
!
!     A      [kg/s]
!     T      [K]
!     b      [kg K/s] 
!     Flux   [kg/s]
!     CT*,   [kg K/s] 
!     DT*,   [kg K/s] 
!     XT*,   [kg K/s]
! 
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  a    => sol % a
  b    => sol % b % val

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! User function
  call User_Mod_Beginning_Of_Compute_Energy(flow, dt, ini)

  ! Old values (o and oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      t % oo(c) = t % o(c)
      t % o (c) = t % n(c)
    end do
  end if

  ! Gradients
  call Grad_Mod_Variable(t, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Energy(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Energy(blend)

  ! Compute tmax and tmin
  if(adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, t % n, t_min, t_max)
  end if

  ! New values
  do c = 1, grid % n_cells
    t % a(c) = 0.0
    t % c(c) = 0.0  ! use t % c for upwind advective fluxes
  end do

  !----------------------------------!
  !   Browse through all the faces   !
  !----------------------------------!
  do s=1,grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ts =      grid % f(s)  * t % n(c1)   &
       + (1.0-grid % f(s)) * t % n(c2)

    ! Compute ts with desired advection scheme
    if(adv_scheme .ne. CENTRAL) then
      call Advection_Scheme(flow, ts, s, t % n, t_min, t_max,  &
                            t % x, t % y, t % z,               &
                            grid % dx, grid % dy, grid % dz,   &
                            adv_scheme, blend) 
    end if

    ! Compute advection term
    if(c2 > 0) then
      t % a(c1) = t % a(c1)-flux(s)*ts*capacity
      t % a(c2) = t % a(c2)+flux(s)*ts*capacity
    else
      t % a(c1) = t % a(c1)-flux(s)*ts*capacity
    end if

    ! Store upwinded part of the advection term in "c"
    if(flux(s).lt.0) then   ! from c2 to c1
      t % c(c1) = t % c(c1)-flux(s)*t % n(c2) * capacity
      if(c2 > 0) then
        t % c(c2) = t % c(c2)+flux(s)*t % n(c2) * capacity
      end if
    else
      t % c(c1) = t % c(c1)-flux(s)*t % n(c1) * capacity
      if(c2 > 0) then
        t % c(c2) = t % c(c2)+flux(s)*t % n(c1) * capacity
      end if
    end if

  end do  ! through faces

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + t % a(c) - t % c(c)
  end do

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  ! Set t % c back to zero 
  do c = 1, grid % n_cells
    t % c(c) = 0.0  
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  if(turbulence_model .ne. NONE .and.  &
     turbulence_model .ne. DNS) then
    call Control_Mod_Turbulent_Prandtl_Number(pr_t)  ! get default pr_t (0.9)
  end if

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(turbulence_model .ne. LES_SMAGORINSKY .and.  &
       turbulence_model .ne. LES_DYNAMIC     .and.  &
       turbulence_model .ne. LES_WALE        .and.  &
       turbulence_model .ne. NONE            .and.  &
       turbulence_model .ne. DNS) then
      pr_t1 = Turbulent_Prandtl_Number(grid, c1)
      pr_t2 = Turbulent_Prandtl_Number(grid, c2)
      pr_t  = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2
    end if

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f1 = grid % fw(s) * t % x(c1) + (1.0-grid % fw(s)) * t % x(c2) 
    ty_f1 = grid % fw(s) * t % y(c1) + (1.0-grid % fw(s)) * t % y(c2)
    tz_f1 = grid % fw(s) * t % z(c1) + (1.0-grid % fw(s)) * t % z(c2)
    tx_f2 = tx_f1
    ty_f2 = ty_f1
    tz_f2 = tz_f1
    if(turbulence_model .ne. NONE .and.  &
       turbulence_model .ne. DNS) then
      con_eff1 =      grid % fw(s)  * (conductivity+capacity*vis_t(c1)/pr_t)  &
               + (1.0-grid % fw(s)) * (conductivity+capacity*vis_t(c2)/pr_t)
      con_t    =      grid % fw(s)  * capacity*vis_t(c1)/pr_t &
               + (1.0-grid % fw(s)) * capacity*vis_t(c2)/pr_t
    else
      con_eff1 = conductivity
    end if

    con_eff2 = con_eff1

    if(turbulence_model .eq. K_EPS        .or.  &
       turbulence_model .eq. K_EPS_ZETA_F .or.  &
       turbulence_model .eq. HYBRID_LES_RANS) then
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALLFL) then
          con_eff1 = con_wall(c1)
          con_eff2 = con_eff1
        end if
      end if
    end if

    ! Total (exact) diffusion flux
    f_ex1 = con_eff1 * (  tx_f1 * grid % sx(s)   &
                        + ty_f1 * grid % sy(s)   &
                        + tz_f1 * grid % sz(s))
    f_ex2 = con_eff2 * (  tx_f2 * grid % sx(s)   &
                        + ty_f2 * grid % sy(s)   &
                        + tz_f2 * grid % sz(s))

    ! Implicit diffusion flux
    f_im1 = con_eff1 * a % fc(s)           &
          * (  tx_f1 * grid % dx(s)      &
             + ty_f1 * grid % dy(s)      &
             + tz_f1 * grid % dz(s) )
    f_im2 = con_eff2 * a % fc(s)           &
          * (  tx_f2 * grid % dx(s)      &
             + ty_f2 * grid % dy(s)      &
             + tz_f2 * grid % dz(s) )

    ! Cross diffusion part
    t % c(c1) = t % c(c1) + f_ex1 - f_im1
    if(c2 > 0) then
      t % c(c2) = t % c(c2) - f_ex2 + f_im2
    end if

    !---------------------------!
    !                           !
    !   Turbulent heat fluxes   !
    !                           !
    !---------------------------!
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
       turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

      ! Turbulent heat fluxes according to GGDH scheme
      ! (first line is GGDH, second line is SGDH substratced 
      ut_s =  (    grid % fw(s)  * ut % n(c1)   &
           +  (1.0-grid % fw(s)) * ut % n(c2))
      vt_s =  (    grid % fw(s)  * vt % n(c1)   &
           +  (1.0-grid % fw(s)) * vt % n(c2))
      wt_s =  (    grid % fw(s)  * wt % n(c1)   &
           +  (1.0-grid % fw(s)) * wt % n(c2))
      t_stress = - (  ut_s * grid % sx(s)                  &
                    + vt_s * grid % sy(s)                  &
                    + wt_s * grid % sz(s) )                &
                    - (con_t * (  tx_f1 * grid % sx(s)     &
                                + ty_f1 * grid % sy(s)     &
                                + tz_f1 * grid % sz(s)) )

      ! Put the influence of turbulent heat fluxes explicitly in the system
      b(c1) = b(c1) + t_stress
      if(c2 > 0) then
        b(c2) = b(c2) - t_stress
      end if
    end if  ! if models are of RSM type

    ! Calculate the coefficients for the sysytem matrix
    a12 = con_eff1 * a % fc(s)
    a21 = con_eff2 * a % fc(s)

    a12 = a12  - min(flux(s), 0.0) * capacity
    a21 = a21  + max(flux(s), 0.0) * capacity

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % dia(c1))  = a % val(a % dia(c1)) + a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) + a21
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
    else if(c2.lt.0) then
      ! Outflow is included because of the flux 
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cell_Type(t, c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cell_Type(t, c2) .eq. CONVECT) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * t % n(c2)
      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALLFL) then
        b(c1) = b(c1) + grid % s(s) * t % q(c2)
      end if
    end if

  end do  ! through sides

  ! Cross diffusion terms are treated explicity
  do c = 1, grid % n_cells
    if(t % c(c) >= 0) then
      b(c)  = b(c) + t % c(c)
    else
      a % val(a % dia(c)) = a % val(a % dia(c))  &
                          - t % c(c) / (t % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Control_Mod_Time_Integration_Scheme(td_scheme)

  ! Two time levels; Linear interpolation
  if(td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * t % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * t % o(c) - 0.5 * a0 * t % oo(c)
    end do
  end if

  call User_Mod_Source(flow, t, a, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for t   !
  !                                 !
  !---------------------------------!

  ! Set under-relaxation factor, then overwrite with control file if specified
  urf = 0.7
  call Control_Mod_Simple_Underrelaxation_For_Energy(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf) * t % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Energy_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of iterations then overwrite with control file if specified
  niter =  5
  call Control_Mod_Max_Iterations_For_Energy_Solver(niter)

  call Bicg(sol,      &
            t % n,    &
            b,        &
            precond,  &
            niter,    &
            tol,      &
            ini_res,  &
            t % res)

  call Info_Mod_Iter_Fill_At(2, 4, t % name, niter, t % res)

  call Comm_Mod_Exchange_Real(grid, t % n)

  ! User function
  call User_Mod_End_Of_Compute_Energy(flow, dt, ini)

  end subroutine
