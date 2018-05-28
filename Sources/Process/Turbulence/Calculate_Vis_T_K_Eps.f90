!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   For k-eps model :                                                          !
!                                                                              !
!   vis_t = c_mu * rho * k^2 * eps                                             !
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   vis_tw = y^+ * vis_t kappa / (E * ln(y^+))                                 !
!                                                                              !
!   For k-eps-v2f model :                                                      !
!                                                                              !
!   vis_t = CmuD * rho * Tsc  * vv                                             !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Control_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Work_Mod, only: re_t => r_cell_01,  &
                      f_mu => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, c
  real    :: pr
  real    :: kin_visc ![m^2/s]
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear    [1/s]        !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t    [kg/(m*s)]   !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity[kg/(m*s)]   !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n  [m^2/s^2]    !
!   Cell volume   vol      [m^3]       | Length          lf       [m]          !
!   left hand s.  A        [kg/s]      | right hand s.   b        [kg*m^2/s^3] !
!   Wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  kin_visc = viscosity/density

  !---------------------!
  !   k-epsilon model   !
  !---------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)
  re_t    = 0.
  f_mu    = 0.
  y_plus  = 0.

  do c = 1, grid % n_cells
    vis_t(c) = c_mu * density * kin % n(c)**2 / (eps % n(c) + TINY)
  end do

  ! Low-Re varaint
  if(turbulence_model_variant .eq. LOW_RE) then
    do c = 1, grid % n_cells
      re_t(c) = kin % n(c)**2. / (kin_visc * eps % n(c) + TINY)
      f_mu(c) = exp(-3.4/(1.0 + 0.02*re_t(c))**2)
      vis_t(c) = f_mu(c) * vis_t(c)
    end do

  ! High-Re varaint
  else
    if(rough_walls .eq. NO) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

            y_plus(c1) = sqrt(tau_wall(c1)/density) * &
              grid % wall_dist(c1) / kin_visc
            vis_wall(c1) = y_plus(c1) * viscosity * kappa / &
              log(e_log * y_plus(c1))
          end if
        end if
      end do
    else if(rough_walls .eq. YES) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALLFL) then
            y_plus(c1) = sqrt(tau_wall(c1)/density) * &
              (grid % wall_dist(c1) + Zo) / kin_visc
            vis_wall(c1) = min(y_plus(c1) * viscosity * kappa  &
              / log((grid % wall_dist(c1) + Zo) / Zo), 1.0e+6 * viscosity)
          end if
        end if
      end do
    end if
  end if

  ! Effective condctivity
  if(heat_transfer .eq. YES) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          pr = viscosity * capacity / conductivity
          con_wall(c1) = viscosity * capacity / pr
        end if
      end if
    end do
  end if

  call Comm_Mod_Exchange(grid, vis_t)

  end subroutine
