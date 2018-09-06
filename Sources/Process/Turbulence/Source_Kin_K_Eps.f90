!==============================================================================!
  subroutine Source_Kin_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
  use Work_Mod, only: kin_x => r_cell_01,  &
                      kin_y => r_cell_02,  &
                      kin_z => r_cell_03,  &
                      kin_sq=> r_cell_04
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
  real    :: u_tot2, u_nor, u_nor2, u_tan
  real    :: kin_vis ![m^2/s]
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  do c = 1, grid % n_cells
    kin_sq(c) = sqrt(kin % n(c))
  end do

  kin_vis = viscosity / density

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  if(turbulence_wall_treatment .eq. HIGH_RE) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

          ! Compute tangential velocity component
          u_tot2 = u % n(c1) **2  &
                 + v % n(c1) **2  &
                 + w % n(c1) **2
          u_nor = ( u % n(c1) * grid % sx(s)     &
                  + v % n(c1) * grid % sy(s)     &
                  + w % n(c1) * grid % sz(s) )   &
                       / sqrt(  grid % sx(s)**2  &
                              + grid % sy(s)**2  &
                              + grid % sz(s)**2 )
          u_nor2 = u_nor*u_nor

          if(u_tot2  >  u_nor2) then
            u_tan = sqrt(u_tot2 - u_nor2)
          else
            u_tan = TINY
          end if

          ! Compute nondimensional wall distance and wall-shear stress
          if(.not. rough_walls) then
            y_plus(c1) = c_mu25 * kin_sq(c1) * grid % wall_dist(c1)  &
                    / kin_vis
            tau_wall(c1) = abs(c_mu25*kappa*density* kin_sq(c1) * u_tan  &
                         / (log(e_log*y_plus(c1))))

            ! Production:
            p_kin(c1) = c_mu25*tau_wall(c1)/density * kin_sq(c1)  &
                      / (kappa*grid % wall_dist(c1))

          else if(rough_walls) then
            y_plus(c1) = c_mu25*kin_sq(c1) * (grid % wall_dist(c1)+Zo)  &
             / kin_vis
            tau_wall(c1) = abs(c_mu25*kappa*density * kin_sq(c1) * u_tan  &
                         / (log((grid % wall_dist(c1)+Zo)/Zo)))

            p_kin(c1) = tau_wall(c1)/density * c_mu25 * kin_sq(c1)  &
                      / (kappa*(grid % wall_dist(c1)+Zo))
            kin % n(c2) = (tau_wall(c1)/density)/0.09**0.5
          end if

          ! Filling up the source term
          b(c1) = b(c1) + density * p_kin(c1) * grid % vol(c1)
        end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
      end if    ! c2 < 0
    end do

    !-----------------------------------------!
    !   Compute the sources in the interior   !
    !-----------------------------------------!
    do c = 1, grid % n_cells
      ! grid % cell_near_wall ensures not to put p_kin twice into &
      ! the near wall cells
      if(.not. grid % cell_near_wall(c)) then
        ! Production:
        p_kin(c)= vis_t(c)/density * shear(c)**2
        b(c) = b(c) + density * p_kin(c) * grid % vol(c)
      end if

      ! Dissipation:
      a % val(a % dia(c)) = a % val(a % dia(c)) + &
          density * eps % n(c)/(kin % n(c)+TINY)*grid % vol(c)
    end do
  end if ! end if mode = wf


  !--------------------------------------------------------!
  !   Jones-Launder model and Launder-Sharma + Yap model   !
  !--------------------------------------------------------!
  if(turbulence_wall_treatment .eq. LOW_RE) then

    do c = 1, grid % n_cells
      p_kin(c) = vis_t(c)/density * shear(c)**2
      ! Production:
      b(c) = b(c) + density * p_kin(c) * grid % vol(c)

      ! Dissipation:
      a % val(a % dia(c)) = a % val(a % dia(c)) + &
        density * eps % n(c)/(kin % n(c)+TINY)*grid % vol(c)
    end do

    call Grad_Mod_For_Phi(grid, kin_sq, 1, kin_x, .true.)  ! dk/dx
    call Grad_Mod_For_Phi(grid, kin_sq, 2, kin_y, .true.)  ! dk/dy
    call Grad_Mod_For_Phi(grid, kin_sq, 3, kin_z, .true.)  ! dk/dz

    do c = 1, grid % n_cells

      a % val(a % dia(c)) = a % val(a % dia(c))                &
                          + 2.0 * viscosity*(  kin_x(c)**2     &
                                             + kin_y(c)**2     &
                                             + kin_z(c)**2 )   &
                           * grid % vol(c) / (kin % n(c) + TINY)
    end do
  end if

  call Comm_Mod_Exchange_Real(grid, kin % n)

  end subroutine
