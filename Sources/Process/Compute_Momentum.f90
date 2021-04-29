!==============================================================================!
  subroutine Compute_Momentum(flow, turb, mult, sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt
  integer, intent(in)           :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Matrix_Type), pointer :: m
  type(Var_Type),    pointer :: ui, uj, uk, t, p
  type(Face_Type),   pointer :: v_flux
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: ui_i(:), ui_j(:), ui_k(:), uj_i(:), uk_i(:)
  real, contiguous,  pointer :: si(:), sj(:), sk(:), di(:), dj(:), dk(:)
  real, contiguous,  pointer :: fi(:), hi(:), cell_fi(:)
  integer                    :: s, c, c1, c2, i
  real                       :: f_ex, f_im, f_stress
  real                       :: vel_max, dt
  real                       :: m0, m12, m21
  real                       :: vis_eff
  real                       :: ui_i_f, ui_j_f, ui_k_f, uj_i_f, uk_i_f
  real                       :: grav_i, p_drop_i
  real                       :: ui_si, ui_di
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
!     [m]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     m              [kg/s]
!     u, v, w        [m/s]
!     bu, bv, bw     [kgm/s^2]   [N]
!     p, pp          [kg/ms^2]   [N/m^2]
!     v_flux         [m^3/s]
!     au*, av*, aw*  [kgm/s^2]   [N]
!     du*, dv*, dw*  [kgm/s^2]   [N]
!     cu*, cv*, cw*  [kgm/s^2]   [N]
!     Wall visc.      vis_w [kg/(m*s)]
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Momentum (without solvers)')

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  t      => flow % t
  p      => flow % p
  dt     =  flow % dt
  m      => sol % m
  b      => sol % b % val

  ! User function
  call User_Mod_Beginning_Of_Compute_Momentum(flow, turb, mult, sol,  &
                                              curr_dt, ini)

  ! Buoyancy force if no VOF is used (VOF calls it from its main)
  ! Now this is a bit stupid.  I do understand that VOF changes only once
  ! before iterations in a time step begin, but buoyancy forces might still
  ! be changing due to temperture differrences which are evoling at the
  ! same time as momentum equations.
  !
  ! On another note, don't call this for the case of VOF because it somehow
  ! fiddles around faces, and that is something which was already done when
  ! computing phyisical properties, so don't have to do it twice.
  if(mult % model .ne. VOLUME_OF_FLUID) then
    call Field_Mod_Buoyancy_Forces(flow)
  end if

  ! Store velocities and volue fluxes from previous iteration
  ! needed for Majumdar's improvement of Rhie and Chow method
  flow % u_star(:)      = flow % u % n(:)
  flow % v_star(:)      = flow % v % n(:)
  flow % w_star(:)      = flow % w % n(:)
  flow % v_flux_star(:) = flow % v_flux % n(:)

  !--------------------------------------------!
  !                                            !
  !                                            !
  !   Browse through all velocity components   !
  !                                            !
  !                                            !
  !--------------------------------------------!
  do i = 1, 3

    if(i .eq. 1) then
      ui   => flow % u;   uj   => flow % v;   uk   => flow % w
      ui_i => ui % x;     ui_j => ui % y;     ui_k => ui % z
      si   => grid % sx;  sj   => grid % sy;  sk   => grid % sz
      di   => grid % dx;  dj   => grid % dy;  dk   => grid % dz
      hi   => p % x;      uj_i => uj % x;     uk_i => uk % x
      fi       => flow % fx
      cell_fi  => flow % cell_fx
      grav_i   =  grav_x
      p_drop_i =  bulk % p_drop_x
    end if
    if(i .eq. 2) then
      ui   => flow % v;   uj   => flow % w;   uk   => flow % u
      ui_i => ui % y;     ui_j => ui % z;     ui_k => ui % x
      si   => grid % sy;  sj   => grid % sz;  sk   => grid % sx
      di   => grid % dy;  dj   => grid % dz;  dk   => grid % dx
      hi   => p % y;      uj_i => uj % y;     uk_i => uk % y
      fi       => flow % fy
      cell_fi  => flow % cell_fy
      grav_i   =  grav_y
      p_drop_i =  bulk % p_drop_y
    end if
    if(i .eq. 3) then
      ui   => flow % w;   uj   => flow % u;   uk   => flow % v
      ui_i => ui % z;     ui_j => ui % x;     ui_k => ui % y
      si   => grid % sz;  sj   => grid % sx;  sk   => grid % sy
      di   => grid % dz;  dj   => grid % dx;  dk   => grid % dy
      hi   => p % z;      uj_i => uj % z;     uk_i => uk % z
      fi       => flow % fz
      cell_fi  => flow % cell_fz
      grav_i   =  grav_z
      p_drop_i =  bulk % p_drop_z
    end if

    ! Initialize matrix and right hand side
    fi     (:) = 0.0  ! all "internal" forces acting on this component
    m % val(:) = 0.0
    b      (:) = 0.0
    fi     (:) = 0.0
    f_stress   = 0.0

    ! Calculate velocity magnitude for normalization
    vel_max = MICRO
    do c = -grid % n_bnd_cells, grid % n_cells
      vel_max = max(vel_max, sqrt(ui % n(c)**2 + uj % n(c)**2 + uk % n(c)**2))
    end do
    call Comm_Mod_Global_Max_Real(vel_max)

    ! Old values (o) and older than old (oo)
    if (flow % piso_status .eqv. .false.) then
      if(ini .eq. 1) then
        do c = 1, grid % n_cells
          ui % oo(c) = ui % o(c)
          ui % o (c) = ui % n(c)
        end do
      end if
    end if

    !---------------!
    !               !
    !   Advection   !
    !               !
    !---------------!
    call Numerics_Mod_Advection_Term(ui, flow % density, v_flux % n, fi)

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

      call Turb_Mod_Face_Vis   (turb, vis_eff,  s)
      call Turb_Mod_Face_Stress(turb, ui, f_stress, s)

      ui_i_f = grid % fw(s)*ui_i(c1) + (1.0-grid % fw(s))*ui_i(c2)
      ui_j_f = grid % fw(s)*ui_j(c1) + (1.0-grid % fw(s))*ui_j(c2)
      ui_k_f = grid % fw(s)*ui_k(c1) + (1.0-grid % fw(s))*ui_k(c2)
      uj_i_f = grid % fw(s)*uj_i(c1) + (1.0-grid % fw(s))*uj_i(c2)
      uk_i_f = grid % fw(s)*uk_i(c1) + (1.0-grid % fw(s))*uk_i(c2)

      ui_si = (  (ui_i_f + ui_i_f) * si(s)    &
               + (ui_j_f + uj_i_f) * sj(s)    &
               + (ui_k_f + uk_i_f) * sk(s) )
      ui_di = (  ui_i_f * di(s)  &
               + ui_j_f * dj(s)  &
               + ui_k_f * dk(s))

      ! Total (exact) viscous stress
      f_ex = vis_eff * ui_si

      ! Implicit viscous stress
      m0 = vis_eff * m % fc(s)
      f_im = ui_di * m0

      ! Cross diffusion part
      ui % c(c1) = ui % c(c1) + f_ex - f_im + f_stress * flow % density(c1)
      if(c2  > 0) then
        ui % c(c2) = ui % c(c2) - f_ex + f_im - f_stress * flow % density(c2)
      end if

      ! Compute the coefficients for the sysytem matrix
      m12 = m0 - min(v_flux % n(s), 0.0) * flow % density(c1)
      m21 = m0 + max(v_flux % n(s), 0.0) * flow % density(c2)

      ! Fill the system matrix
      if(c2 > 0) then
        m % val(m % pos(1,s)) = m % val(m % pos(1,s)) - m12
        m % val(m % dia(c1))  = m % val(m % dia(c1))  + m12
        m % val(m % pos(2,s)) = m % val(m % pos(2,s)) - m21
        m % val(m % dia(c2))  = m % val(m % dia(c2))  + m21
      else if(c2  < 0) then
        ! Outflow is not included because it was causing problems
        if((Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)  .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)    .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL)) then
           ! (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) ) then
          m % val(m % dia(c1)) = m % val(m % dia(c1)) + m12
          fi(c1) = fi(c1) + m12 * ui % n(c2)
        end if
      end if

      ! Here we clean up momentum from the false diffusion
      call Turb_Mod_Substract_Face_Stress(turb, ui_si, ui_di,            &
                                                ui % n(c1), ui % n(c2),  &
                                                m % fc(s), fi, s)

    end do  ! through faces

    ! Explicit treatment for cross diffusion terms
    do c = 1, grid % n_cells
      fi(c) = fi(c) + ui % c(c)
    end do

    !--------------------!
    !                    !
    !   Inertial terms   !
    !                    !
    !--------------------!
    call Numerics_Mod_Inertial_Term(ui, flow % density, m, fi, dt)

    !---------------------------------!
    !                                 !
    !   Various force contributions   !
    !                                 !
    !---------------------------------!

    !--------------------------!
    !   Global pressure drop   !
    !--------------------------!
    do c = 1, grid % n_cells
      fi(c) = fi(c) + p_drop_i * grid % vol(c)
    end do

    !--------------------!
    !   Buoyancy force   !  (Units are ugly, different from others :-/)
    !--------------------!
    do c = 1, grid % n_cells
      fi(c) = fi(c) + cell_fi(c)
    end do

    !----------------------------------!
    !   Surface tension contribution   !
    !----------------------------------!
    if(mult % model .eq. VOLUME_OF_FLUID) then
      call Multiphase_Mod_Vof_Momentum_Contribution(mult, fi, i)
    end if

    !----------------------------------------!
    !   All other terms defined by the user  !
    !----------------------------------------!
    call User_Mod_Force(flow, ui, m, fi)

    !-----------------------------------------------------------!
    !                                                           !
    !   Copy forces from current component to right hand side   !
    !   (Note: pressure gradients are not with other forces)    !
    !                                                           !
    !-----------------------------------------------------------!
    do c = 1, grid % n_cells
      b(c) = fi(c) - hi(c) * grid % vol(c)
    end do

    !------------------------------------------------!
    !   Save the coefficients from the discretized   !
    !   momentum equation before under-relaxation    !
    !------------------------------------------------!
    do c = 1, grid % n_cells
      m % sav(c) = m % val(m % dia(c))
    end do

    !----------------------------------------------!
    !   Explicit solution for the PISO algorithm   !
    !----------------------------------------------!
    call Compute_Momentum_Explicit(flow, ui, sol)

    !-----------------------------------!
    !                                   !
    !   Solve the equations for u,v,w   !
    !                                   !
    !-----------------------------------!

    !--------------------------------------------------------!
    !   If not inside the PRIME part of the PISO algorithm   !
    !--------------------------------------------------------!
    if(flow % piso_status .eqv. .false.) then

      ! Under-relax the equations
      call Numerics_Mod_Under_Relax(ui, m, b)

      ! Call linear solver
      call Cpu_Timer_Mod_Start('Linear_Solver_For_Momentum')

      call Solver_Mod_Bicg(sol,           &
                           m,             &
                           ui % n,        &
                           b,             &
                           ui % precond,  &
                           ui % mniter,   &
                           ui % eniter,   &
                           ui % tol,      &
                           ui % res,      &
                           norm = vel_max)
      call Cpu_Timer_Mod_Stop('Linear_Solver_For_Momentum')

      ! Fill the info screen up
      if (flow % p_m_coupling == SIMPLE) then
        call Info_Mod_Iter_Fill_At(1, i, ui % name, ui % eniter, ui % res)
      end if

    end if

  end do  ! browsing through components

  !------------------------------------------------------------------------!
  !   Save the coefficients from the discretized momentum equation after   !
  !   under-relaxation.  This is the way it should be done, but analysis   !
  !   (analythical and numerical) shows that it only changes the conver-   !
  !   gence history.  If saving was engaged here, matrix for pressure      !
  !   would have smaller entries, pressure correction would be higher,     !
  !   which would, in theory, lead to higher velocity corrections and      !
  !   faster convergence.  But, since many cases in T-Flows are already    !
  !   tuned for coefficients saved befoe under-relaxation, I am hesitant   !
  !   to abruptly make a change in the code and estimate them here.        !
  !------------------------------------------------------------------------!
  !@ do c = 1, grid % n_cells
  !@   m % sav(c) = m % val(m % dia(c))
  !@ end do

  ! Refresh buffers for m % sav before discretizing for pressure
  call Grid_Mod_Exchange_Cells_Real(grid, m % sav)

  ! User function
  call User_Mod_End_Of_Compute_Momentum(flow, turb, mult, sol, curr_dt, ini)

  call Cpu_Timer_Mod_Stop('Compute_Momentum (without solvers)')

  end subroutine
