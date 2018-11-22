!==============================================================================!
  subroutine Source_Eps_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in the eps transport equation,                   !
!   wall shear stress (wall function approuch)                                 !
!------------------------------------------------------------------------------!
!   int( density (c_1e eps/kin Gk - c_2e eps^2/kin) )dV                        !
!                                                                              !
!   assigns epsilon from the wall function:                                    !
!                                                                              !
!   Eps_w = Cmu^(3/4)* Kin^(3/2)/(Kappa/y)                                     !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Roughness_Coefficient 
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, j
  real    :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real    :: re_t, f_mu
  real    :: eps_wf, eps_int, ebf, y_star, u_tau_new, fa, kin_vis
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  kin_vis = viscosity/density

  do c = 1, grid % n_cells
    ! Positive contribution:
    b(c) = b(c) + &
            c_1e * p_kin(c) * eps % n(c)/kin % n(c) * grid % vol(c)

    ! Negative contribution:
    re_t = kin % n(c)*kin % n(c)/(kin_vis*eps % n(c))
    y_star = sqrt(sqrt(kin_vis * eps % n(c))) *     &
             grid % wall_dist(c)/kin_vis
    f_mu = (1.0 - exp(-y_star/3.1))**2              &
         * (1.0 - 0.3*exp(-(re_t/6.5)*(re_t/6.5)))

    f_mu = min(f_mu,1.0)

    a % val(a % dia(c)) = a % val(a % dia(c)) &
     + density * f_mu* c_2e * eps % n(c) / kin % n(c) * grid % vol(c)
  end do

  ! Imposing a boundary condition on wall for eps
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
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
        u_nor_sq = u_nor**2

        if( u_tot_sq  > u_nor_sq) then
          u_tan = sqrt(u_tot_sq - u_nor_sq)
        else
          u_tan = TINY
        end if
 
        if(rough_walls) then 
          z_o = Roughness_Coefficient(grid, c1)       
          eps % n(c1) = c_mu75 * kin % n(c1)**1.5  &
                      / ((grid % wall_dist(c1) + z_o) * kappa)

          ! Adjusting coefficient to fix eps value in near wall calls
          do j = a % row(c1), a % row(c1 + 1) - 1
            a % val(j) = 0.0 
          end do

          b(c1) = eps % n(c1) * density
          a % val(a % dia(c1)) = 1.0 * density
        else
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan   &
                       / log(e_log * max(y_plus(c1), 1.05))

          u_tau_new = sqrt(tau_wall(c1)/density)
          y_plus(c1) = Y_Plus_Low_Re(u_tau_new, grid % wall_dist(c1), kin_vis)
          ebf = 0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1))

          eps_int = 2.0*viscosity/density * kin % n(c1)    &
                  / grid % wall_dist(c1)**2
          eps_wf  = c_mu75 * kin % n(c1)**1.5              &
                  / (grid % wall_dist(c1) * kappa)

          if(y_plus(c1) > 3) then
            fa = min(density*u_tau_new**3  &
               / (kappa*grid % wall_dist(c1)*p_kin(c1)),1.0)

            eps % n(c1) = (1.0-fa)*eps_int + fa*eps_wf

            ! Adjusting coefficient to fix eps value in near wall calls
            do j = a % row(c1), a % row(c1 + 1) - 1
              a % val(j) = 0.0
            end do

            b(c1) = eps % n(c1)
            a % val(a % dia(c1)) = 1.0
          else
            eps % n(c2) = eps_int
          end if ! y_plus > 4
        end if   ! rough_walls
      end if     ! wall or wall_flux
    end if       ! c2 < 0
  end do

  end subroutine
