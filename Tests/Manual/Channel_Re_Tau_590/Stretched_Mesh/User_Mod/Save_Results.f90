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
  type(Var_Type),  pointer :: u, v, w, kin, eps, zeta, f22, omega, t
  integer                  :: n_prob, i, c, count, n_points, fu, c1, c2, s
  character(SL)            :: res_name
  real, allocatable        :: z_p(:), u_p(:), uw_p(:),                    &
                              kin_p(:), eps_p(:), f22_p(:),               &
                              zeta_p(:), omg_p(:),                        &
                              y_plus_p(:), t_p(:), ind(:), wall_p(:)
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p, t_tau
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const, t_wall
  logical                  :: there
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Flow % Alias_Energy  (t)   
  omega => Turb % omega

  ! Read constant physical properties from control file
  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  ubulk    = bulk % flux_x / bulk % area_x
  n_points = 0

  ! Set file names to save results                                               
  call File % Set_Name(res_name,                      &                          
                       time_step = Time % Curr_Dt(),  &                          
                       extension = '-res.dat')                                   
                                                                                 
  ! Number of probes in cell rows                                                
  n_prob = Grid % n_z_planes - 1                  

  !-------------------------------------------------------------------------!
  !   Allocate memory for variables to be extracted in homogeneous planes   !
  !-------------------------------------------------------------------------!
  allocate(n_p     (n_prob));  n_p      = 0
  allocate(wall_p  (n_prob));  wall_p   = 0.0
  allocate(u_p     (n_prob));  u_p      = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(t_p (n_prob));      t_p      = 0.0
  allocate(f22_p   (n_prob));  f22_p    = 0.0
  allocate(omg_p   (n_prob));  omg_p    = 0.0
  allocate(zeta_p  (n_prob));  zeta_p   = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob  ! n_prob is in cell rows                                     
    do c = Cells_In_Domain()                                                     
      if(Grid % zc(c) > (Grid % z_coord_plane(i)) .and.  &                       
         Grid % zc(c) < (Grid % z_coord_plane(i+1))) then   

        wall_p  (i) = wall_p  (i) + Grid % wall_dist(c)
        u_p     (i) = u_p     (i) + u % n(c)
        uw_p    (i) = uw_p    (i) + Turb % vis_t(c) * (u % z(c) + w % x(c))
        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)
        if(Turb % model == K_EPS_ZETA_F) then
          kin_p   (i) = kin_p   (i) + kin % n(c)
          eps_p   (i) = eps_p   (i) + eps % n(c)
          f22_p   (i) = f22_p   (i) + f22  % n(c)
          zeta_p  (i) = zeta_p  (i) + zeta % n(c)
        else if(Turb % model == K_EPS) then
          kin_p   (i) = kin_p   (i) + kin % n(c)
          eps_p   (i) = eps_p   (i) + eps % n(c)
        else if(Turb % model == K_OMEGA_SST) then
          kin_p   (i) = kin_p   (i) + kin % n(c)
          omg_p   (i) = omg_p   (i) + omega % n(c)
        end if 
        if(Flow % heat_transfer) t_p (i) = t_p (i) + t  % n(c)     

        n_count(i) = n_count(i) + 1
      end if
    end do
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      wall_p  (i) = wall_p  (i) / n_count(i)
      u_p     (i) = u_p     (i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      omg_p   (i) = omg_p   (i) / n_count(i)
      uw_p    (i) = uw_p    (i) / n_count(i)
      t_p (i)     = t_p (i)     / n_count(i)
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

  if(Flow % heat_transfer) then
    t_wall = 0.0
    t_tau     = Flow % heat_flux / (dens_const * capa_const * u_tau_p)
    do s = 1, Grid % n_faces                                                     
      c1 = Grid % faces_c(1,s)                                                   
      c2 = Grid % faces_c(2,s)                                                   
      if(c2  < 0) then                                                           
        if( Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &                           
            Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then                           
                                                                                 
          t_wall  = t_wall + t % n(c2)                                           
          n_points = n_points + 1                                                
        end if                                                                   
      end if                                                                     
    end do                                                                       
                                                                                 
    call Global % Sum_Real(t_wall)                                               
    call Global % Sum_Int(n_points)                                              
                                                                                 
    call Global % Wait                                                           
                                                                                 
    t_wall  = t_wall / n_points    
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


  do i = 1, n_prob
    wall_p(i) = dens_const * wall_p(i) * u_tau_p / visc_const
    u_p   (i) = u_p  (i) / u_tau_p
    kin_p (i) = kin_p(i) / u_tau_p**2
    eps_p (i) = eps_p(i)*visc_const / (u_tau_p**4*dens_const)
    uw_p  (i) = uw_p (i) / (u_tau_p**2 * dens_const)
    f22_p (i) = f22_p(i) * visc_const / u_tau_p**2
    if(Flow % heat_transfer) t_p (i) = (t_wall - t_p(i)) / t_tau
  end do

  if(Turb % model == K_EPS_ZETA_F) then
    write(fu,'(a)') '#  1:z,  2:u,  3:kin,  4:eps,  5:uw,'    //  &
                    '  6:f22,  7:zeta,  8:vis_t/visc_const'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(8es15.5e3)')  wall_p(i),   &  !  1
                                 u_p(i),      &  !  2
                                 kin_p(i),    &  !  3
                                 eps_p(i),    &  !  4
                                 uw_p(i),     &  !  5
                                 f22_p(i),    &  !  6
                                 zeta_p(i),   &  !  7
                                 t_p(i)          !  8
      end if
    end do
  else if(Turb % model == K_EPS) then
    write(fu,'(a)') '#  1:z,  2:u,  3:kin,  4:eps,  5:uw,'    //  &
                    '   6:vis_t/visc_const'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(6es15.5e3)')  wall_p(i),   &  !  1
                                 u_p(i),      &  !  2
                                 kin_p(i),    &  !  3
                                 eps_p(i),    &  !  4
                                 uw_p(i),     &  !  5
                                 t_p(i)          !  6
      end if
    end do
  else if(Turb % model == K_OMEGA_SST) then
    write(fu,'(a)') '#  1:z,  2:u,  3:kin,  4:omega,  5:uw,'    //  &
                    '   6:vis_t/visc_const'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(6es15.5e3)')  wall_p(i),   &  !  1
                                 u_p(i),      &  !  2
                                 kin_p(i),    &  !  3
                                 omg_p(i),    &  !  4
                                 uw_p(i),     &  !  5
                                 t_p(i)          !  6
      end if
    end do
  else 
    write(fu,'(a)') '#  1:z,  2:u' 
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(3es15.5e3)')  wall_p(i),   &  !  1
                                 u_p(i),      &  !  2
                                 t_p(i)          !  3
      end if
    end do
  end if
  
  close(fu)

  if(First_Proc())  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
