!==============================================================================!
  subroutine Time_And_Length_Scale(Turb, Grid)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Grid_Type),  target :: Grid
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  real                      :: kin_vis   ! kinematic viscosity [m^2/s]
  integer                   :: c
  real, contiguous, pointer :: t_1(:), t_2(:), t_3(:), l_1(:), l_2(:), l_3(:)
  real, contiguous, pointer :: eps_l(:)
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          l_1       [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!==============================================================================!

  call Work % Connect_Real_Cell(t_1, t_2, t_3, l_1, l_2, l_3, eps_l)

  ! Take aliases
  Flow => Turb % pnt_flow
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)

  if(Turb % model .eq. K_EPS_ZETA_F) then

    do c = Cells_In_Domain_And_Buffers()
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      kin_vis = Flow % viscosity(c) / Flow % density(c)

      t_1(c) = kin % n(c) / eps_l(c)
      t_2(c) = c_t*sqrt(kin_vis/eps_l(c))
      t_3(c) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(c) * Flow % shear(c) + TINY)

      l_1(c) = kin % n(c)**1.5 / eps_l(c)
      l_2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25
      l_3(c) = sqrt(kin % n(c)/3.0)  &
             / (c_mu_d * zeta % n(c) * Flow % shear(c) + TINY)

      Turb % t_scale(c) =       max( min(t_1(c), t_3(c)), t_2(c) )
      Turb % l_scale(c) = c_l * max( min(l_1(c), l_3(c)), l_2(c) )
    end do

  else if(Turb % model .eq. HYBRID_LES_RANS) then
    do c = Cells_In_Domain_And_Buffers()
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      kin_vis = Flow % viscosity(c) / Flow % density(c)

      t_1(c) = kin % n(c) / eps_l(c)
      t_2(c) = c_t*sqrt(kin_vis/eps_l(c))
      t_3(c) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(c) * Flow % shear(c) + TINY)

      l_1(c) = kin % n(c)**1.5 / eps_l(c)
      l_2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25

      Turb % t_scale(c) =       max( min(t_1(c), t_3(c)), t_2(c) )
      Turb % l_scale(c) = c_l * max( l_1(c), l_2(c) )
    end do

  else if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then

    do c = Cells_In_Domain_And_Buffers()
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      kin_vis = Flow % viscosity(c) / Flow % density(c)

      t_1(c) = kin % n(c)/eps_l(c)
      t_2(c) = c_t*sqrt(kin_vis/eps_l(c))

      l_1(c) = kin % n(c)**1.5/eps_l(c)
      l_2(c) = c_nu * (kin_vis**3 / eps_l(c))**0.25

      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)

      Turb % t_scale(c) =       max( t_1(c), t_2(c) )
      Turb % l_scale(c) = c_l * max( l_1(c), l_2(c) )
    end do

  else if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    do c = Cells_In_Domain_And_Buffers()
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      t_1(c) = kin % n(c)/eps_l(c)
      l_1(c) = kin % n(c)**1.5/eps_l(c)

      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
      Turb % t_scale(c) = t_1(c)
      Turb % l_scale(c) = l_1(c)
    end do

  else if(Turb % model .eq. K_EPS) then

    do c = Cells_In_Domain_And_Buffers()
      eps_l(c) = eps % n(c) + TINY ! limited eps % n

      t_1(c) = kin % n(c) / eps_l(c)
      t_3(c) = 0.6/(sqrt(6.0) * c_mu * Flow % shear(c) + TINY)
      Turb % t_scale(c) =  min(t_1(c), t_3(c))
    end do
  end if

  call Work % Disconnect_Real_Cell(t_1, t_2, t_3, l_1, l_2, l_3, eps_l)

  end subroutine
