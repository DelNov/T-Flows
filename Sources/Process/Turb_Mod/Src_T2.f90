!==============================================================================!
  subroutine Turb_Mod_Src_T2(turb, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in t2 transport equation for k-eps_t2 model      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, eps, ut, vt, wt
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum  (flow, u, v, w)
  call Field_Mod_Alias_Energy    (flow, t)
  call Turb_Mod_Alias_K_Eps      (turb, kin, eps)
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Solver_Mod_Alias_System   (sol, a, b)

  call Grad_Mod_Array(grid, t % n, t % x, t % y, t % z, .true.)

  !-----------------------------------------!
  !   Compute the sources in all the cells  !
  !-----------------------------------------!

  ! Production source:
  do c = 1, grid % n_cells

    turb % p_t2(c) = - 2.0 * (  ut % n(c) * t % x(c)   &
                              + vt % n(c) * t % y(c)   &
                              + wt % n(c) * t % z(c))

    b(c) = b(c) + turb % p_t2(c) * grid % vol(c)

   ! Negative contribution
   a % val(a % dia(c)) = a % val(a % dia(c)) +  &
         2.0 * density(c)  * eps % n(c) / (kin % n(c) + TINY) * grid % vol(c)

  end do

  end subroutine

