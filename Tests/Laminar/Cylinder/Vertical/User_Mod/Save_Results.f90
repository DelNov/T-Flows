!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  integer, optional         :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Var_Type),  pointer :: kin, eps, zeta, f22
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: ut, vt, wt
  integer                  :: s, c1, c2, n_points
  real                     :: nuss_mean, t_wall, alfa_const, ra, pr, d
  real                     :: nuss_cc, h_cc, q_cc
  real                     :: rho_const, mu_const, nu_const
  real                     :: capa_const, k_const
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)

  ! Take constant physical properties
  call Control % Mass_Density        (rho_const)
  call Control % Dynamic_Viscosity   (mu_const)
  nu_const = mu_const / rho_const
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(k_const)

  d      =  1.0  ! characteristic dim.
  t_wall = 50.0  ! temp. at the wal

  alfa_const = k_const / (capa_const * rho_const)  ! thermal diffusivity
  ra = (t_wall - Flow % t_ref) * d**3  &
     / (nu_const * alfa_const)                     ! rayleigh number
  pr = nu_const / alfa_const                       ! prandtl number

  ! Churchil and Chu formula
  h_cc = k_const * (0.6 + 0.387 * ra ** (1.0/6.0)                         &
                        / (1.0 + (0.559/pr) ** (9.0/16.0)) ** (8.0/27.0)  &
                      ) ** 2.0
  q_cc = h_cc * Flow % heated_area * (t_wall - Flow % t_ref)

  ! Nusselt number from Churchil and Chu
  nuss_cc = h_cc * d / k_const

  nuss_mean = 0.0
  n_points  = 0

  if(Flow % heat_transfer) then
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2  < 0) then
        if( Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
            Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

          nuss_mean = nuss_mean + t % q(c2)  &
                    / (k_const * (t % n(c2) - Flow % t_ref + TINY))
          n_points = n_points + 1
        end if
      end if
    end do

    call Global % Sum_Real(nuss_mean)
    call Global % Sum_Int(n_points)

    call Global % Wait

    nuss_mean = nuss_mean / n_points
  end if

  if(First_Proc()) then
    print *, 't_ref       = ', Flow % t_ref
    print *, 'alfa_const  = ', alfa_const
    print *, 'ra          = ', ra
    print *, 'pr          = ', pr
    print *, 'h_cc        = ', h_cc
    print *, 'q_cc        = ', q_cc
    print *, 'nuss_cc     = ', nuss_cc
    print *, 'heat        = ', Flow % heat
    print *, 'heat_flux   = ', Flow % heat_flux
    print *, 'heated_area = ', Flow % heated_area
    print *, 'nuss_mean   = ', nuss_mean
  end if

  end subroutine
