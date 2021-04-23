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
                      u_f => r_face_01,  &
                      v_f => r_face_02,  &
                      w_f => r_face_03
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
  real,              pointer :: u_relax
  real                       :: a12, px_f, py_f, pz_f, dens_h, fs
  integer                    :: c, s, c1, c2
!==============================================================================!

  call Cpu_Timer_Mod_Start('Rhie_And_Chow')

  ! Take aliases
  grid    => flow % pnt_grid
  p       => flow % p
  v_flux  => flow % v_flux
  u_relax => flow % u_rel_corr
  a       => sol % a
  m       => sol % m
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  ! call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, ini)

  !------------------------------!
  !   Work out cell velocities   !
  !------------------------------!
  do c = 1, grid % n_cells
    u_c(c) = flow % u % n(c)
    v_c(c) = flow % v % n(c)
    w_c(c) = flow % w % n(c)
  end do

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
      u_f(s) = u_c(c1) * fs + u_c(c2) * (1.0-fs)
      v_f(s) = v_c(c1) * fs + v_c(c2) * (1.0-fs)
      w_f(s) = w_c(c1) * fs + w_c(c2) * (1.0-fs)

      ! Calculate coeficients for the system matrix
      ! a12 [m*m^3*s/kg = m^4s/kg]
      a12 = u_relax * 0.5 * a % fc(s)                    &
                    * ( grid % vol(c1) / m % sav(c1)     &
                      + grid % vol(c2) / m % sav(c2) )

      ! Interpolate pressure gradients as proposed by Denner
      ! (Equation 3.57 in his PhD thesis)
      ! dens_h           [kg/m^3]
      ! px_f, py_f, pz_f [kg/m^2s^2]
      dens_h = 2.0 / (1.0 / flow % density(c1) + 1.0 / flow % density(c2))
      px_f = 0.5 * dens_h * (  p % x(c1) / flow % density(c1)  &
                             + p % x(c2) / flow % density(c2) )
      py_f = 0.5 * dens_h * (  p % y(c1) / flow % density(c1)  &
                             + p % y(c2) / flow % density(c2) )
      pz_f = 0.5 * dens_h * (  p % z(c1) / flow % density(c1)  &
                             + p % z(c2) / flow % density(c2) )

      ! Calculate current volume flux through cell face with pressure
      ! defined at a cell face and assuming that pressure correction
      ! (pp) part is treated implicitly
      v_flux % n(s) = u_f(s) * grid % sx(s)          &
                    + v_f(s) * grid % sy(s)          &
                    + w_f(s) * grid % sz(s)          &
                    + a12 * (p % n(c1) - p % n(c2))  &
                    + a12 * (  px_f * grid % dx(s)   &
                             + py_f * grid % dy(s)   &
                             + pz_f * grid % dz(s))

    end if

  end do

  call Cpu_Timer_Mod_Stop('Rhie_And_Chow')

  end subroutine
