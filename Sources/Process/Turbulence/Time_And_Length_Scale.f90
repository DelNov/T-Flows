!==============================================================================!
  subroutine Time_And_Length_Scale(grid)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
  use Work_Mod, only: t1    => r_cell_12,  &  ! [s]
                      t2    => r_cell_13,  &  ! [s]
                      t3    => r_cell_14,  &  ! [s]
                      l1    => r_cell_15,  &  ! [m]
                      l2    => r_cell_16,  &  ! [m]
                      l3    => r_cell_17,  &  ! [m]
                      eps_l => r_cell_18      ! [m]
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  real    :: kin_vis   ! kinematic viscosity [m^2/s]
  integer :: c
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!

  kin_vis = viscosity / density

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    do c = 1, grid % n_cells
      eps_l(c) = max(eps % n(c),tiny) 
 
      t1(c) = kin % n(c)/eps_l(c)
      t2(c) = c_t*sqrt(kin_vis/eps_l(c))
      t3(c) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(c) * shear(c) + TINY)

      l1(c) = kin % n(c)**1.5/eps_l(c)
      l2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25
      l3(c) = sqrt(kin % n(c)/3.0)/(c_mu_d * zeta % n(c) * shear(c) + TINY)
     
      t_scale(c) =       max( min(t1(c), t3(c)), t2(c) )
      l_scale(c) = c_l * max( min(l1(c), l3(c)), l2(c) )
    end do

  else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

    do c = 1, grid % n_cells
      eps_l(c) = max(eps % n(c),tiny) 
      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
 
      t1(c) = kin % n(c)/eps_l(c)
      t2(c) = c_t*sqrt(kin_vis/eps_l(c))

      l1(c) = kin % n(c)**1.5/eps_l(c)
      l2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25

      t_scale(c) =       max( t1(c), t2(c) )
      l_scale(c) = c_l * max( l1(c), l2(c) )
    end do

  else if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    do c = 1, grid % n_cells
      eps_l(c) = max(eps % n(c),tiny) 
      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
 
      t1(c) = kin % n(c)/eps_l(c)

      l1(c) = kin % n(c)**1.5/eps_l(c)

      t_scale(c) = t1(c)
      l_scale(c) = l1(c)
    end do

  end if

  end subroutine
