!==============================================================================!
  subroutine Time_And_Length_Scale(grid, turb)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Field_Mod
  use Turb_Mod
  use Grid_Mod
  use Var_Mod
  use Work_Mod, only: t_1   => r_cell_12,  &  ! [s]
                      t_2   => r_cell_13,  &  ! [s]
                      t_3   => r_cell_14,  &  ! [s]
                      l_1   => r_cell_15,  &  ! [m]
                      l_2   => r_cell_16,  &  ! [m]
                      l_3   => r_cell_17,  &  ! [m]
                      eps_l => r_cell_18      ! [m]
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type), target :: grid
  type(Turb_Type), target :: turb
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  real                      :: kin_vis   ! kinematic viscosity [m^2/s]
  integer                   :: c
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          l_1       [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    do c = 1, grid % n_cells
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      kin_vis = flow % viscosity(c) / flow % density(c)

      t_1(c) = kin % n(c) / eps_l(c)
      t_2(c) = c_t*sqrt(kin_vis/eps_l(c))
      t_3(c) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(c) * flow % shear(c) + TINY)

      l_1(c) = kin % n(c)**1.5 / eps_l(c)
      l_2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25
      l_3(c) = sqrt(kin % n(c)/3.0)  &
             / (c_mu_d * zeta % n(c) * flow % shear(c) + TINY)

      turb % t_scale(c) =       max( min(t_1(c), t_3(c)), t_2(c) )
      turb % l_scale(c) = c_l * max( min(l_1(c), l_3(c)), l_2(c) )
    end do

  else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

    do c = 1, grid % n_cells
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      kin_vis = flow % viscosity(c) / flow % density(c)

      t_1(c) = kin % n(c)/eps_l(c)
      t_2(c) = c_t*sqrt(kin_vis/eps_l(c))

      l_1(c) = kin % n(c)**1.5/eps_l(c)
      l_2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25

      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)

      turb % t_scale(c) =       max( t_1(c), t_2(c) )
      turb % l_scale(c) = c_l * max( l_1(c), l_2(c) )
    end do

  else if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    do c = 1, grid % n_cells
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      t_1(c) = kin % n(c)/eps_l(c)
      l_1(c) = kin % n(c)**1.5/eps_l(c)

      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
      turb % t_scale(c) = t_1(c)
      turb % l_scale(c) = l_1(c)
    end do

  end if

  end subroutine
