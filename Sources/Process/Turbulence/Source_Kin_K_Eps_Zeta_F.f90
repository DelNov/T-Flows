!==============================================================================!
  subroutine Source_Kin_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation.                       !
!------------------------------------------------------------------------------!
!   In kinetic energy eq. there are two source terms:                          !
!   int( density (p_kin - eps ) )dV                                            !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  integer :: c, c1, c2, s
  real    :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real    :: lf
  real    :: alpha1, l_rans, l_sgs
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Production source:
  do c = 1, grid % n_cells
    p_kin(c) = vis_t(c)/density * shear(c)**2.
    b(c)     = b(c) + density * p_kin(c) * grid % vol(c)
  end do

  if(turbulence_model .eq. K_EPS_ZETA_F .and.  &
     turbulence_statistics) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      l_sgs  = 0.8*lf
      l_rans = 0.41*grid % wall_dist(c)
      alpha1 = max(1.0,l_rans/l_sgs)

      if(alpha1 < 1.05) then
        a % val(a % dia(c)) = a % val(a % dia(c)) + &
             density * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)
      else
        a % val(a % dia(c)) = a % val(a % dia(c)) +  &
          density *  min(alpha1**1.45 * eps % n(c), kin % n(c)**1.5  &
          / (lf*0.01)) / (kin % n(c) + TINY) * grid % vol(c)
      end if
    end do
  else
    do c = 1, grid % n_cells

      a % val(a % dia(c)) = a % val(a % dia(c)) + &
           density * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)

      if(buoyancy) then
        buoy_beta(c) = 1.0
        g_buoy(c) = -buoy_beta(c) * (grav_x * ut % n(c) +  &
                                     grav_y * vt % n(c) +  &
                                     grav_z * wt % n(c)) * density
        b(c) = b(c) + max(0.0, g_buoy(c) * grid % vol(c))
        a % val(a % dia(c)) = a % val(a % dia(c))  &
                  + max(0.0,-g_buoy(c) * grid % vol(c) / (kin % n(c) + TINY))
      end if
    end do
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tot_sq = u % n(c1) * u % n(c1) &
                 + v % n(c1) * v % n(c1) &
                 + w % n(c1) * w % n(c1)
        u_nor  = ( u % n(c1) * grid % sx(s)     &
                 + v % n(c1) * grid % sy(s)     &
                 + w % n(c1) * grid % sz(s) )   &
                 / sqrt(  grid % sx(s)*grid % sx(s)  &
                        + grid % sy(s)*grid % sy(s)  &
                        + grid % sz(s)*grid % sz(s))
        u_nor_sq = u_nor*u_nor

        if( u_tot_sq  > u_nor_sq) then
          u_tan = sqrt(u_tot_sq - u_nor_sq)
        else
          u_tan = TINY
        end if

        if(y_plus(c1) > 3.0) then
          if(.not. rough_walls) then
            ! Wall shear s.
            tau_wall(c1) = density*kappa*u_tau(c1)*u_tan / &
              (log(e_log*y_plus(c1)))
            p_kin(c1) = tau_wall(c1)/density * u_tau(c1) / &
              (kappa*grid % wall_dist(c1))
          else if (rough_walls) then
            tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                           /(log((grid % wall_dist(c1)+Zo)/Zo))
            p_kin(c1) = tau_wall(c1)/density * u_tau(c1) / &
                        (kappa*(grid % wall_dist(c1)+Zo))
            kin % n(c2) = (tau_wall(c1)/density) / 0.09**0.5
          end if

          b(c1) = b(c1) + &
            ( density*p_kin(c1) - vis_t(c1)*shear(c1)**2 ) * grid % vol(c1)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  end subroutine
