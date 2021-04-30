!==============================================================================!
  subroutine Rhie_And_Chow(flow, mult, sol)
!------------------------------------------------------------------------------!
!   Computes face velocitites with Rhie and Chow interpolation method          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: u_c => r_cell_01,  &
                      v_c => r_cell_02,  &
                      w_c => r_cell_03,  &
                      v_m => r_cell_04,  &  ! for Rhie and Chow
                      t_m => r_cell_05,  &  ! for Choi
                      u_f => r_face_01,  &
                      v_f => r_face_02,  &
                      w_f => r_face_03
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence must be established                  !
!                                                                              !
!   Main_Pro                        (allocates Work_Mod)                       !
!     |                                                                        |
!     +----> Compute_Pressure       (doesn't use Work_Mod)                     !
!             |                                                                |
!             +----> Rhie_And_Chow  (safe to use Work_Mod)                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: a               ! pressure matrix
  type(Matrix_Type), pointer :: m               ! momentum matrix
  real                       :: a12, px_f, py_f, pz_f, fs, dens_h
  integer                    :: s, c1, c2, c
  logical, parameter         :: CHOI = .false.
!==============================================================================!

  call Cpu_Timer_Mod_Start('Rhie_And_Chow')

  ! Take aliases
  grid   => flow % pnt_grid
  p      => flow % p
  v_flux => flow % v_flux
  a      => sol % a
  m      => sol % m
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  ! call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, ini)

  !--------------------------------------!
  !   Take velocities as last computed   !
  !--------------------------------------!
  do c = 1, grid % n_cells
    u_c(c) = flow % u % n(c)
    v_c(c) = flow % v % n(c)
    w_c(c) = flow % w % n(c)
  end do

  !--------------------------------------!
  !   Store grid % vol(c) / m % sav(c)   !
  !--------------------------------------!
  ! Units here: m^3 s / kg
  do c = 1, grid % n_cells
    v_m(c) = grid % vol(c) / m % sav(c)
  end do

  ! First part (cell-centered) of Choi's correction
  ! (Subtract the cell-centered unsteady terms)
  if(CHOI) then
    do c = 1, grid % n_cells

      ! Unit for t_m: m^3 * kg/m^3 / s * s/kg = 1
      t_m(c) = (grid % vol(c) * flow % density(c) / flow % dt) / m % sav(c)

      u_c(c) = u_c(c) - t_m(c) * u % o(c)
      v_c(c) = v_c(c) - t_m(c) * v % o(c)
      w_c(c) = w_c(c) - t_m(c) * w % o(c)
    end do
  end if

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if(c2 > 0) then

      ! Interpolate velocity 
      u_f(s) = fs * u_c(c1) + (1.0 - fs) * u_c(c2)
      v_f(s) = fs * v_c(c1) + (1.0 - fs) * v_c(c2)
      w_f(s) = fs * w_c(c1) + (1.0 - fs) * w_c(c2)

      ! Calculate coeficients for the pressure matrix
      ! Units: m * m^3 s / kg = m^4 s / kg
      a12 = 0.5 * a % fc(s) * (v_m(c1) + v_m(c2))
      a % val(a % pos(1,s)) = -a12
      a % val(a % pos(2,s)) = -a12
      a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) +  a12

      ! Interpolate pressure gradients
      ! Units: kg/(m^2 s^2) * m^3 s / kg * m = m^2 / s
      px_f = 0.5 * (p % x(c1)*v_m(c1) + p % x(c2)*v_m(c2)) * grid % dx(s)
      py_f = 0.5 * (p % y(c1)*v_m(c1) + p % y(c2)*v_m(c2)) * grid % dy(s)
      pz_f = 0.5 * (p % z(c1)*v_m(c1) + p % z(c2)*v_m(c2)) * grid % dz(s)

      ! Calculate mass or volume flux through cell face
      ! Units in lines which follow:
      ! m^3 / s = m/s * m^2
      !         + m^4 s / kg * kg / (m s^2)
      !         + m * m^2/s = m^3/s
      v_flux % n(s) = (  u_f(s) * grid % sx(s)             &
                       + v_f(s) * grid % sy(s)             &
                       + w_f(s) * grid % sz(s) )           &
                       + a12 * (p % n(c1) - p % n(c2))     &
                       + a % fc(s) * (px_f + py_f + pz_f)

      ! Second part of Choi's correction
      ! (Add face-centered flux from previous time step)
      if(CHOI) then
        v_flux % n(s) = v_flux % n(s)  &
                      + v_flux % o(s) * 0.5 * (t_m(c1) + t_m(c2))
      end if

    end if
  end do

  call Cpu_Timer_Mod_Stop('Rhie_And_Chow')

  end subroutine
