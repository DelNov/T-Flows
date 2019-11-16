!==============================================================================!
  subroutine Compute_Momentum(flow, turb, mult, i, sol, dt, ini)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Cpu_Timer_Mod,  only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use Field_Mod,      only: Field_Type, buoyancy, t_ref,  &
                            grav_x, grav_y, grav_z,       &
                            density, viscosity
  use Turb_Mod
  use Multiphase_Mod, only: Multiphase_Type, &
                            Multiphase_Mod_Vof_Surface_Tension_Contribution,  &
                            multiphase_model, VOLUME_OF_FLUID
  use Var_Mod,        only: Var_Type
  use Face_Mod,       only: Face_Type
  use Grid_Mod,       only: Grid_Type
  use Bulk_Mod,       only: Bulk_Type
  use Info_Mod,       only: Info_Mod_Iter_Fill_At
  use Numerics_Mod
  use Solver_Mod,     only: Solver_Type, Solver_Mod_Alias_System, Bicg, Cg, Cgs
  use Matrix_Mod,     only: Matrix_Type
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  integer                       :: i           ! component
  type(Solver_Type),     target :: sol
  real                          :: dt
  integer                       :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Matrix_Type), pointer :: a
  type(Var_Type),    pointer :: ui, uj, uk, t, p
  type(Face_Type),   pointer :: m_flux
  real,              pointer :: b(:)
  real,              pointer :: ui_i(:), ui_j(:), ui_k(:), uj_i(:), uk_i(:)
  real,              pointer :: si(:), sj(:), sk(:), di(:), dj(:), dk(:)
  real,              pointer :: h_i(:)
  integer                    :: s, c, c1, c2, exec_iter
  real                       :: f_ex, f_im, f_stress
  real                       :: vel_max
  real                       :: a0, a12, a21
  real                       :: vis_eff, vis_ts
  real                       :: ui_i_f, ui_j_f, ui_k_f, uj_i_f, uk_i_f
!------------------------------------------------------------------------------!
!
!  Stress tensor on the face s:
!
!    T = mi * [    2*dU/dx     dU/dy+dV/dx   dU/dz+dW/dx  ]
!             [  dU/dy+dV/dx     2*dV/dy     dV/dz+dW/dy  ]
!             [  dU/dz+dW/dx   dV/dz+dW/dy     2*dW/dz    ]
!
!  The forces, acting on the cell face are:
!
!    Fx = T11*Sx + T12*Sy + T13*Sz
!    Fy = T21*Sx + T22*Sy + T23*Sz
!    Fz = T31*Sx + T32*Sy + T33*Sz
!
!  which could also be written in the compact form:
!
!    {F} = [T]{S}
!
!  U component:
!
!    Fx = Txx*Sx + Txy*Sy + Txz*Sz
!
!  V component:
!
!    Fy = Tyx*Sx + Tyy*Sy + Tyz*Sz
!
!  W component:
!
!    Fz = Tzx*Sx + Tzy*Sy + Tzz*Sz
!
!------------------------------------------------------------------------------!
!
!  The form of equations which I am solving:
!
!     /             /              /               /
!    |     du      |              |               |
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | p dS
!    |     dt      |              |               |
!   /             /              /               /
!
!  Dimension of the system under consideration
!
!     [a]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     a              [kg/s]
!     u, v, w        [m/s]
!     bu, bv, bw     [kgm/s^2]   [N]
!     p, pp          [kg/ms^2]   [N/m^2]
!     flux           [kg/s]
!     au*, av*, aw*  [kgm/s^2]   [N]
!     du*, dv*, dw*  [kgm/s^2]   [N]
!     cu*, cv*, cw*  [kgm/s^2]   [N]
!     Wall visc.      vis_w [kg/(m*s)]
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Momentum (without solvers)')

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  m_flux => flow % m_flux
  t      => flow % t
  p      => flow % p
  call Solver_Mod_Alias_System(sol, a, b)

  if(i .eq. 1) then
    ui   => flow % u;   uj   => flow % v;   uk   => flow % w
    ui_i => ui % x;     ui_j => ui % y;     ui_k => ui % z
    si   => grid % sx;  sj   => grid % sy;  sk   => grid % sz
    di   => grid % dx;  dj   => grid % dy;  dk   => grid % dz
    h_i  => p % x;      uj_i => uj % x;     uk_i => uk % x
  end if
  if(i .eq. 2) then
    ui   => flow % v;   uj   => flow % w;   uk   => flow % u
    ui_i => ui % y;     ui_j => ui % z;     ui_k => ui % x
    si   => grid % sy;  sj   => grid % sz;  sk   => grid % sx
    di   => grid % dy;  dj   => grid % dz;  dk   => grid % dx
    h_i  => p % y;      uj_i => uj % y;     uk_i => uk % y
  end if
  if(i .eq. 3) then
    ui   => flow % w;   uj   => flow % u;   uk   => flow % v
    ui_i => ui % z;     ui_j => ui % x;     ui_k => ui % y
    si   => grid % sz;  sj   => grid % sx;  sk   => grid % sy
    di   => grid % dz;  dj   => grid % dx;  dk   => grid % dy
    h_i  => p % z;      uj_i => uj % z;     uk_i => uk % z
  end if

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0
  f_stress   = 0.0

  ! User function
  call User_Mod_Beginning_Of_Compute_Momentum(flow, dt, ini)

  ! Calculate velocity magnitude for normalization
  vel_max = 0.0
  do c = -grid % n_bnd_cells, grid % n_cells
    vel_max = max(vel_max, sqrt(ui % n(c)**2 + uj % n(c)**2 + uk % n(c)**2))
  end do
  call Comm_Mod_Global_Max_Real(vel_max)

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      ui % oo(c) = ui % o(c)
      ui % o (c) = ui % n(c)
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(ui, 1.0, m_flux % n, sol,  &
                                   ui_i,                      &
                                   ui_j,                      &
                                   ui_k,                      &
                                   di,                        &
                                   dj,                        &
                                   dk)

  !---------------!
  !               !
  !   Diffusion   !
  !               !
  !---------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    call Turb_Mod_Calculate_Face_Vis   (turb, vis_eff,  s)
    call Turb_Mod_Calculate_Face_Stress(turb, ui, f_stress, s)

    ui_i_f = grid % fw(s)*ui_i(c1) + (1.0-grid % fw(s))*ui_i(c2)
    ui_j_f = grid % fw(s)*ui_j(c1) + (1.0-grid % fw(s))*ui_j(c2)
    ui_k_f = grid % fw(s)*ui_k(c1) + (1.0-grid % fw(s))*ui_k(c2)
    uj_i_f = grid % fw(s)*uj_i(c1) + (1.0-grid % fw(s))*uj_i(c2)
    uk_i_f = grid % fw(s)*uk_i(c1) + (1.0-grid % fw(s))*uk_i(c2)

    ! Total (exact) viscous stress
    f_ex = vis_eff*(      2.0*ui_i_f  * si(s)      &
                    + (ui_j_f+uj_i_f) * sj(s)      &
                    + (ui_k_f+uk_i_f) * sk(s) )

    a0 = vis_eff * a % fc(s)

    ! Implicit viscous stress
    f_im = (   ui_i_f*di(s)                &
             + ui_j_f*dj(s)                &
             + ui_k_f*dk(s)) * a0

    ! Cross diffusion part
    ui % c(c1) = ui % c(c1) + f_ex - f_im + f_stress * density(c1)
    if(c2  > 0) then
      ui % c(c2) = ui % c(c2) - f_ex + f_im - f_stress * density(c2)
    end if

    ! Compute the coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(m_flux % n(s), real(0.0))
    a21 = a21  + max(m_flux % n(s), real(0.0))

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a21
    else if(c2  < 0) then
      ! Outflow is not included because it was causing problems
      if((Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)  .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)    .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL)) then
         ! (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1) = b(c1) + a12 * ui % n(c2)
      end if
    end if

  end do  ! through faces

  ! Here we clean up momentum from the false diffusion
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant .ne. STABILIZED) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        vis_ts =     (grid % fw(s)  * turb % vis_t(c1)  &
               + (1.0-grid % fw(s)) * turb % vis_t(c2))
        a0 = a % fc(s)*vis_ts
        vis_eff = vis_ts

        ui_i_f = grid % fw(s) * ui_i(c1) + (1.0-grid % fw(s)) * ui_i(c2)
        ui_j_f = grid % fw(s) * ui_j(c1) + (1.0-grid % fw(s)) * ui_j(c2)
        ui_k_f = grid % fw(s) * ui_k(c1) + (1.0-grid % fw(s)) * ui_k(c2)
        uj_i_f = grid % fw(s) * uj_i(c1) + (1.0-grid % fw(s)) * uj_i(c2)
        uk_i_f = grid % fw(s) * uk_i(c1) + (1.0-grid % fw(s)) * uk_i(c2)

        f_ex = vis_eff*( 2.0*ui_i_f         * si(s) &
                          + (ui_j_f+uj_i_f) * sj(s) &
                          + (ui_k_f+uk_i_f) * sk(s) )

        f_im = (  ui_i_f * di(s)  &
                + ui_j_f * dj(s)  &
                + ui_k_f * dk(s)) * vis_eff * a % fc(s)

        b(c1) = b(c1) - vis_eff * (ui % n(c2) -ui % n(c1)) * a % fc(s)  &
              - f_ex + f_im
        if(c2  > 0) then
          b(c2) = b(c2) + vis_eff * (ui % n(c2) -ui % n(c1)) * a % fc(s)  &
                + f_ex - f_im
        end if
      end do
    end if
  end if

  ! Explicit treatment for cross diffusion terms
  do c = 1, grid % n_cells
    b(c) = b(c) + ui % c(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  call Numerics_Mod_Inertial_Term(ui, density, sol, dt)

  !---------------------------------!
  !                                 !
  !   Pressure term contributions   !
  !                                 !
  !---------------------------------!

  !--------------------------!
  !   Global pressure drop   !
  !--------------------------!
  if(ui % name .eq. 'U') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_x * grid % vol(c)
    end do
  else if(ui % name .eq. 'V') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_y * grid % vol(c)
    end do
  else if(ui % name .eq. 'W') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_z * grid % vol(c)
    end do
  end if

  !---------------------------------!
  !   Local pressure distribution   !
  !---------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) - h_i(c) * grid % vol(c)
  end do

  !--------------------!
  !   Buoyancy force   !
  !--------------------!
  if(buoyancy) then
    if(ui % name .eq. 'U') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density(c) * grav_x * (t % n(c) - t_ref)  &
             * flow % beta * grid % vol(c)
      end do
    else if(ui % name .eq. 'V') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density(c) * grav_y * (t % n(c) - t_ref)  &
             * flow % beta * grid % vol(c)
      end do
    else if(ui % name .eq. 'W') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density(c) * grav_z * (t % n(c) - t_ref)  &
             * flow % beta * grid % vol(c)
      end do
    end if

  !-------------------!
  !   Gravity force   !
  !-------------------!
  else
    if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) > TINY) then
      if(ui % name .eq. 'U') then
        do c = 1, grid % n_cells
          b(c) = b(c) + density(c) * grav_x * grid % vol(c)
        end do
      else if(ui % name .eq. 'V') then
        do c = 1, grid % n_cells
          b(c) = b(c) + density(c) * grav_y * grid % vol(c)
        end do
      else if(ui % name .eq. 'W') then
        do c = 1, grid % n_cells
          b(c) = b(c) + density(c) * grav_z * grid % vol(c)
        end do
      end if
    end if
  end if

  !----------------------------------!
  !   Surface tension contribution   !
  !----------------------------------!
  if(multiphase_model .eq. VOLUME_OF_FLUID) then

    if (surface_tension > TINY ) then
      if(ui % name .eq. 'U') then
        call Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
        call Grad_Mod_Variable(mult % vof)
        do c = 1, grid % n_cells
          b(c) = b(c) + surface_tension * mult % vof % oo(c)  &
                                        * mult % vof % x(c)   &
                                        * grid % vol(c)
        end do
      else if(ui % name .eq. 'V') then
        do c = 1, grid % n_cells
          b(c) = b(c) + surface_tension * mult % vof % oo(c)  &
                                        * mult % vof % y(c)   &
                                        * grid % vol(c)
        end do
      else if(ui % name .eq. 'W') then
        do c = 1, grid % n_cells
          b(c) = b(c) + surface_tension * mult % vof % oo(c)  &
                                        * mult % vof % z(c)   &
                                        * grid % vol(c)
        end do
      end if

    end if

  end if

  !----------------------------------------!
  !   All other terms defined by the user  !
  !----------------------------------------!
  call User_Mod_Force(grid, ui, a, b)

  ! Save the coefficients for pressure equation
  do c = 1, grid % n_cells
    a % sav(c) = a % val(a % dia(c))
  end do

  !-----------------------------------!
  !                                   !
  !   Solve the equations for u,v,w   !
  !                                   !
  !-----------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(ui, sol)

  ! Call linear solver
  call Cpu_Timer_Mod_Start('Linear_Solver_For_Momentum')
  call Bicg(sol,           &
            ui % n,        &
            b,             &
            ui % precond,  &
            ui % niter,    &
            exec_iter,     &
            ui % tol,      &
            ui % res,      &
            norm = vel_max)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Momentum')

  ! Fill the info screen up
  if(ui % name .eq. 'U') then
    call Info_Mod_Iter_Fill_At(1, 1, ui % name, exec_iter, ui % res)
  end if
  if(ui % name .eq. 'V') then
    call Info_Mod_Iter_Fill_At(1, 2, ui % name, exec_iter, ui % res)
  end if
  if(ui % name .eq. 'W') then
    call Info_Mod_Iter_Fill_At(1, 3, ui % name, exec_iter, ui % res)
  end if

  call Comm_Mod_Exchange_Real(grid, ui % n)

  ! User function
  call User_Mod_End_Of_Compute_Momentum(flow, dt, ini)

  call Cpu_Timer_Mod_Stop('Compute_Momentum (without solvers)')

  end subroutine
