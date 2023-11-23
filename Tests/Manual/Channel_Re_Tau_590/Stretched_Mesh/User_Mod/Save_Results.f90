!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in the file name_res_plus.dat                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, optional        :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, kin, eps, zeta, f22
  integer                  :: n_prob, i, c, count, n_points, fu
  character(SL)            :: coord_name, res_name
  real, allocatable        :: z_p(:), u_p(:), uw_p(:),                    &
                              kin_p(:), eps_p(:), f22_p(:), zeta_p(:),    &
                              y_plus_p(:),  vis_t_p(:), ind(:), wall_p(:)
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  logical                  :: there
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  ! This version of the Save_Results works only for sequential runs
  if(Parallel_Run()) then
    if(First_Proc()) then
      print *, '#==========================================================='
      print *, '# This version of User_Mod_Save_Results is not intended for '
      print *, '# parallel runs.  Skipped writing of profiles in .dat file! '
      print *, '#-----------------------------------------------------------'
    end if
    return
  end if

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)

  ! Read constant physical properties from control file
  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  ! Set file name for results
  call File % Set_Name(res_name,                      &
                       time_step = Time % Curr_Dt(),  &
                       appendix='-res-plus',          &
                       extension='.dat')

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file=coord_name, exist=there)
  if(.not. there) then
    if(First_Proc()) then
      print *, '#=============================================================='
      print *, '# In order to extract profiles and write them in ascii files'
      print *, '# the code has to read cell-faces coordinates '
      print *, '# in wall-normal direction in the ascii file ''case_name.1d.'''
      print *, '# The file format should be as follows:'
      print *, '# 10  ! number of cells + 1'
      print *, '# 1 0.0'
      print *, '# 2 0.1'
      print *, '# 3 0.2'
      print *, '# ... '
      print *, '#--------------------------------------------------------------'
    end if

    return
  end if

  ubulk    = bulk % flux_x / bulk % area_x
  n_points = 0

  ! Open ASCII file for writing in the first available unit (fu)
  call File % Open_For_Reading_Ascii(coord_name, fu)

  ! Read number of probes ...
  read(fu, *) n_prob

  ! ... allocate memory for them
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! ... and read the intervals positions
  do i = 1, n_prob
    read(fu, *) ind(i), z_p(i)
  end do
  close(fu)

  !-------------------------------------------------------------------------!
  !   Allocate memory for variables to be extracted in homogeneous planes   !
  !-------------------------------------------------------------------------!
  allocate(n_p     (n_prob));  n_p      = 0
  allocate(wall_p  (n_prob));  wall_p   = 0.0
  allocate(u_p     (n_prob));  u_p      = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(vis_t_p (n_prob));  vis_t_p  = 0.0
  allocate(f22_p   (n_prob));  f22_p    = 0.0
  allocate(zeta_p  (n_prob));  zeta_p   = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p  (i) = wall_p  (i) + Grid % wall_dist(c)
        u_p     (i) = u_p     (i) + u % n(c)
        kin_p   (i) = kin_p   (i) + kin % n(c)
        eps_p   (i) = eps_p   (i) + eps % n(c)
        uw_p    (i) = uw_p    (i) + Turb % vis_t(c) * (u % z(c) + w % x(c))
        vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)
        f22_p   (i) = f22_p   (i) + f22  % n(c)
        zeta_p  (i) = zeta_p  (i) + zeta % n(c)

        n_count(i) = n_count(i) + 1
      end if
    end do
  end do

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p  (i) = wall_p  (i) / n_count(i)
      u_p     (i) = u_p     (i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      uw_p    (i) = uw_p    (i) / n_count(i)
      vis_t_p (i) = vis_t_p (i) / n_count(i)
      f22_p   (i) = f22_p   (i) / n_count(i)
      zeta_p  (i) = zeta_p  (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
    end if
  end do

  !------------------------------------!
  !   Non-dimensionalize the results   !
  !------------------------------------!

  ! Calculating friction velocity
  if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z)) / dens_const)
  else
    u_tau_p =  sqrt( (visc_const*sqrt(u_p(1)**2)         &
                                      / wall_p(1))       &
                                      / dens_const)
  end if

  !----------------------------------------!
  !   Write the results in the .dat file   !
  !----------------------------------------!
  call File % Open_For_Writing_Ascii(res_name, fu)

  pr = visc_const * capa_const / cond_const
  re = dens_const * ubulk * 2.0 / visc_const
  cf_dean = 0.073*(re)**(-0.25)
  cf      = u_tau_p**2/(0.5*ubulk**2)
  error   = abs(cf_dean - cf)/cf_dean * 100.0
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Density  = ', dens_const
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', dens_const * ubulk * 2.0 / visc_const
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', dens_const * u_tau_p / visc_const
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p / ubulk)**2
  write(fu,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p
  write(fu,'(a1,(a12,f12.6,a2,a22))') &
  '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'

  write(fu,'(a)') '#  1:z,  2:u,  3:kin,  4:eps,  5:uw,'    //  &
                  '  6:f22,  7:zeta,  8:vis_t/visc_const'

  do i = 1, n_prob-1
    wall_p(i) = dens_const * wall_p(i) * u_tau_p / visc_const
    u_p   (i) = u_p  (i) / u_tau_p
    kin_p (i) = kin_p(i) / u_tau_p**2
    eps_p (i) = eps_p(i)*visc_const / (u_tau_p**4*dens_const)
    uw_p  (i) = uw_p (i) / (u_tau_p**2 * dens_const)
    f22_p (i) = f22_p(i) * visc_const / u_tau_p**2
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(fu,'(8es15.5e3)')  wall_p(i),   &  !  1
                               u_p(i),      &  !  2
                               kin_p(i),    &  !  3
                               eps_p(i),    &  !  4
                               uw_p(i),     &  !  5
                               f22_p(i),    &  !  6
                               zeta_p(i),   &  !  7
                               vis_t_p(i)      !  8
    end if
  end do

  close(fu)

  if(First_Proc())  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
