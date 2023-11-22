!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine was derived from many simular one spread throughout the    !
!   Tests' subdirectories, but is tailored to provide inlet profiles as the    !
!   boundary conditions for the city benchmark cases.                          !
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
  type(Var_Type),  pointer :: u, v, w
  type(Var_Type),  pointer :: kin, eps, zeta, f22
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  integer                  :: n_prob, pl, c, i, s, c1, c2, n_points
  character(SL)            :: coord_name, res_name
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:), uw_p(:),      &
                              kin_p(:), eps_p(:), f22_p(:), zeta_p(:),      &
                              y_plus_p(:),  vis_t_p(:), ind(:), wall_p(:)
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: t_wall, t_tau, d_wall, nu_mean, t_inf
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  logical                  :: there
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  Assert(Turb % model .eq. K_EPS_ZETA_F)

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)

  ! Take constant physical properties
  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  ! Set file name for results
  call File % Set_Name(res_name,                      &
                       time_step = Time % Curr_Dt(),  &
                       appendix  = '-res',            &
                       extension = '.dat')

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
      print *, '# 9 0.2'
      print *, '# ... '
      print *, '#--------------------------------------------------------------'
    end if

    return
  end if

  ubulk    = bulk % flux_x / (bulk % area_x)
  t_wall   = 0.0
  nu_mean  = 0.0
  n_points = 0

  open(9, file=coord_name)

  ! Write the number of searching intervals
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(9,*) ind(pl), z_p(pl)
  end do
  close(9)

  allocate(n_p     (n_prob));  n_p      = 0
  allocate(wall_p  (n_prob));  wall_p   = 0.0
  allocate(u_p     (n_prob));  u_p      = 0.0
  allocate(v_p     (n_prob));  v_p      = 0.0
  allocate(w_p     (n_prob));  w_p      = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(vis_t_p (n_prob));  vis_t_p  = 0.0
  allocate(f22_p   (n_prob));  f22_p    = 0.0
  allocate(zeta_p  (n_prob));  zeta_p   = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + Grid % wall_dist(c)
        u_p   (i) = u_p(i) + u % n(c)
        v_p   (i) = v_p(i) + v % n(c)
        w_p   (i) = w_p(i) + w % n(c)

        if(Turb % model .ne. NO_TURBULENCE_MODEL) then
          kin_p   (i) = kin_p   (i) + kin % n(c)
          eps_p   (i) = eps_p   (i) + eps % n(c)
          uw_p    (i) = uw_p    (i) + Turb % vis_t(c) * (u % z(c) + w % x(c))
          vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
          y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)
        end if

        if(Turb % model .eq. K_EPS_ZETA_F) then
          f22_p (i) = f22_p (i) + f22  % n(c)
          zeta_p(i) = zeta_p(i) + zeta % n(c)
        end if

        n_count(i) = n_count(i) + 1
      end if
    end do
  end do

  ! Average over all processors
  do pl=1, n_prob-1
    call Global % Sum_Int(n_count(pl))

    call Global % Sum_Real(wall_p(pl))

    call Global % Sum_Real(u_p(pl))
    call Global % Sum_Real(v_p(pl))
    call Global % Sum_Real(w_p(pl))

    call Global % Sum_Real(kin_p   (pl))
    call Global % Sum_Real(eps_p   (pl))
    call Global % Sum_Real(uw_p    (pl))
    call Global % Sum_Real(vis_t_p (pl))
    call Global % Sum_Real(y_plus_p(pl))

    call Global % Sum_Real(f22_p (pl))
    call Global % Sum_Real(zeta_p(pl))

  end do

  call Global % Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p  (i) = wall_p(i) / n_count(i)
      u_p     (i) = u_p   (i) / n_count(i)
      v_p     (i) = v_p   (i) / n_count(i)
      w_p     (i) = w_p   (i) / n_count(i)

      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      uw_p    (i) = uw_p    (i) / n_count(i)
      vis_t_p (i) = vis_t_p (i) / n_count(i)
      f22_p   (i) = f22_p   (i) / n_count(i)
      zeta_p  (i) = zeta_p  (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
    end if
  end do

  ! Calculating friction velocity
  if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z)) / dens_const)
  else
    u_tau_p =  sqrt( (visc_const*sqrt(u_p(1)**2 +        &
                                      v_p(1)**2 +        &
                                      w_p(1)**2)         &
                                      / wall_p(1))       &
                                      / dens_const)
  end if

  if(u_tau_p .eq. 0.0) then
    if(First_Proc()) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if

    return
  end if

  open(9, file = res_name)

  pr = visc_const * capa_const / cond_const
  re = dens_const * ubulk * 2.0 / visc_const
  cf = u_tau_p**2/(0.5*ubulk**2)
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'density  = ', dens_const 
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'U        = ', ubulk 
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'Re       = ', dens_const * ubulk * 2.0 / visc_const
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'Re_tau   = ', dens_const * u_tau_p / visc_const
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p / ubulk)**2
  write(9,'(a1,(a12,es15.5e3))')  &
  '#', 'U_tau    = ', u_tau_p

  write(9, '(a19)') '# Number of points:'
  write(9, '(i6)')   n_prob
  write(9, '(a86)') '# The following line can be coppied to control file ' // &
                    ' (without the leading # of course)'
  write(9, '(a80)') '# z            '   //  &
                    '  u            '   //  &
                    '  kin          '   //  &
                    '  eps          '   //  &
                    '  zeta         '   //  &
                    '  f22'

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(9,'(8es15.5e3)')  wall_p(i),   &  !  1
                              u_p(i),      &  !  2
                              kin_p(i),    &  !  3
                              eps_p(i),    &  !  4
                              zeta_p(i),   &  !  5
                              f22_p(i)        !  6
    end if
  end do

  close(9)

  if(First_Proc())  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
