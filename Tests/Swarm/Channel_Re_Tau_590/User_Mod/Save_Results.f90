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
  type(Field_Type),      target :: Flow
  type(Turb_Type),       target :: Turb
  type(Vof_Type),        target :: Vof
  type(Swarm_Type),      target :: Swarm
  integer, optional, intent(in) :: domain    ! current domain
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: n_prob, pl, c, i
  integer                  :: fu1, fu2
  character(SL)            :: coord_name, result_name, result_name_plus
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:),          &
                              ind(:),  wall_p(:), kin_p(:), eps_p(:),  &
                              uw_p(:), uu_p(:), vv_p(:), ww_p(:),      &
                              zeta_p(:), f22_p(:), uw_mod_p(:),        &
                              ww_mod_p(:), y_plus_p(:), vis_t_p(:)
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: ubulk, err, re, cf_dean, cf, u_tau_p
  real                     :: visc_const, dens_const
  logical                  :: there
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  ! Take constant physical properties
  call Control % Mass_Density     (dens_const)
  call Control % Dynamic_Viscosity(visc_const)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

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

  do c = 1, Grid % n_cells
    ubulk    = bulk % flux_x / bulk % area_x
  end do 

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
  allocate(uu_p    (n_prob));  uu_p     = 0.0
  allocate(vv_p    (n_prob));  vv_p     = 0.0
  allocate(ww_p    (n_prob));  ww_p     = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(zeta_p  (n_prob));  zeta_p   = 0.0
  allocate(f22_p   (n_prob));  f22_p    = 0.0
  allocate(uw_mod_p(n_prob));  uw_mod_p = 0.0
  allocate(ww_mod_p(n_prob));  ww_mod_p = 0.0
!  allocate(ww_mod_p(Grid % n_cells));  ww_mod_p = 0.0
  allocate(vis_t_p (n_prob));  vis_t_p  = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0

  !!=========================================================
  !! DEBUGGING 
  !do c = 1, Grid % n_cells - Grid % comm % n_buff_cells
  !  ww_mod_p(c) =  Turb % kin_mean(c) * Turb % zeta_mean(c)
  !end do 
  !if(First_Proc()) then 
  ! print *, "w_mod(503) = ", ww_mod_p(503)
  ! print *, "min_val_wmod = ", minval(ww_mod_p(:))
  ! print *, "max_val_wmod = ", maxval(ww_mod_p(:))
  ! stop  
  !end if
  !!=========================================================


  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, Grid % n_cells - Grid % comm % n_buff_cells 
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + Grid % wall_dist(c)
        u_p   (i) = u_p   (i) + Turb % u_mean(c)
        v_p   (i) = v_p   (i) + Turb % v_mean(c)
        w_p   (i) = w_p   (i) + Turb % w_mean(c)

        uu_p(i) = uu_p(i) + Turb % uu_res(c)  &
                          - Turb % u_mean(c) * Turb % u_mean(c)
        vv_p(i) = vv_p(i) + Turb % vv_res(c)  &
                          - Turb % v_mean(c) * Turb % v_mean(c)
        ww_p(i) = ww_p(i) + Turb % ww_res(c)  &
                          - Turb % w_mean(c) * Turb % w_mean(c)
        uw_p(i) = uw_p(i) + Turb % uw_res(c)  &
                          - Turb % u_mean(c) * Turb % w_mean(c)

        ! Modeled quantities
        kin_p   (i) = kin_p   (i) + Turb % kin_mean(c)
        eps_p   (i) = eps_p   (i) + Turb % eps_mean(c)
        zeta_p  (i) = zeta_p  (i) + Turb % zeta_mean(c)
        f22_p   (i) = f22_p   (i) + Turb % f22_mean(c)
        uw_mod_p(i) = uw_mod_p(i) + Turb % vis_t_eff(c)*(u % y(c) + v % x(c))
        ww_mod_p(i) = ww_mod_p(i) + Turb % kin_mean(c) * Turb % zeta_mean(c)
        vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c) 
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

    call Global % Sum_Real(uu_p(pl))
    call Global % Sum_Real(vv_p(pl))
    call Global % Sum_Real(ww_p(pl))
    call Global % Sum_Real(uw_p(pl))

    call Global % Sum_Real(uw_mod_p(pl))
    call Global % Sum_Real(kin_p   (pl))
    call Global % Sum_Real(eps_p   (pl))
    call Global % Sum_Real(zeta_p  (pl))
    call Global % Sum_Real(f22_p   (pl))
    call Global % Sum_Real(vis_t_p (pl))
    call Global % Sum_Real(y_plus_p(pl))
    call Global % Sum_Real(ww_mod_p(pl))
  end do

  call Global % Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count(i)
      u_p   (i) = u_p   (i) / n_count(i)
      v_p   (i) = v_p   (i) / n_count(i)
      w_p   (i) = w_p   (i) / n_count(i)

      uu_p(i) = uu_p(i) / n_count(i)
      vv_p(i) = vv_p(i) / n_count(i)
      ww_p(i) = ww_p(i) / n_count(i)
      uw_p(i) = uw_p(i) / n_count(i)

      uw_mod_p(i) = uw_mod_p(i) / n_count(i)
      ww_mod_p(i) = ww_mod_p(i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      zeta_p  (i) = zeta_p  (i) / n_count(i)
      f22_p   (i) = f22_p   (i) / n_count(i)
      vis_t_p (i) = vis_t_p (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
    end if
  end do

  ! Calculating friction velocity and friction temperature
!    do c= 1, Grid % n_cells
 if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z))/dens_const)
  else
      u_tau_p = sqrt( (visc_const*sqrt(u_p(1)**2 +   &
                                       v_p(1)**2 +   &
                                       w_p(1)**2)    &
                                       / wall_p(1))  &
                                       / dens_const)
!    end do
  end if

  if(u_tau_p < TINY) then
    if(First_Proc()) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if
    return
  end if

  call File % Set_Name(result_name, time_step = Time % Curr_Dt(),  &
       appendix='-res', extension='.dat')
  call File % Open_For_Writing_Ascii(result_name, fu1)
  call File % Set_Name(result_name_plus, time_step = Time % Curr_Dt(),  &
       appendix='-res-plus', extension='.dat')
  call File % Open_For_Writing_Ascii(result_name_plus, fu2)

  open(fu1,file=result_name)
  open(fu2,file=result_name_plus)

  re = dens_const * ubulk * 2.0 / visc_const
  cf_dean = 0.073*(re)**(-0.25)
  cf      = u_tau_p**2/(0.5*ubulk**2)
  err     = abs(cf_dean - cf)/cf_dean * 100.0

  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk 
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', dens_const * ubulk * 2.0/visc_const
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', dens_const*u_tau_p/visc_const
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
  write(fu1,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p 
  write(fu1,'(a1,(a12,f12.6,a2,a22))') & 
  '#', 'Cf_error = ', err, ' %', 'Dean formula is used.'


  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk 
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', dens_const * ubulk * 2.0/visc_const
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', dens_const*u_tau_p/visc_const
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
  write(fu2,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p 
  write(fu2,'(a1,(a12,f12.6,a2,a22))') & 
  '#', 'Cf_error = ', err, ' %', 'Dean formula is used.'

  write(fu1,'(a1,2x,a105)') '#',  ' 1) z,'                                //  &
                                  ' 2) u mean,'                           //  &
                                  ' 3) kin_resolved, 4) kin_modeled,'     //  &
                                  ' 5) kin_tot,  6) uw_resolved,'         //  &
                                  ' 7) uw_modeled, 8) uw_tot, 9) vis_t,'  //  &
                                  ' 10) ww_total '

  write(fu2,'(a1,2x,a105)') '#',  ' 1) z,'                                //  &
                                  ' 2) u mean,'                           //  &
                                  ' 3) kin_resolved, 4) kin_modeled,'     //  &
                                  ' 5) kin_tot,  6) uw_resolved,'         //  &
                                  ' 7) uw_modeled, 8) uw_tot, 9) vis_t,'  //  &
                                  ' 10) ww_total '

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(fu1,'(10e15.7)')                                         &
                        wall_p(i),                                   &  !  1
                        u_p(i),                                      &  !  2
                        0.5*(uu_p(i)+vv_p(i)+ww_p(i)),               &  !  3
                        kin_p(i),                                    &  !  4
                        (0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i)),  &  !  5
                        uw_p(i),                                     &  !  6
                        uw_mod_p(i),                                 &  !  7
                        (uw_p(i) + uw_mod_p(i)),                     &  !  8
                        vis_t_p(i),                                  &  !  9
                        (ww_p(i) + ww_mod_p(i))                         ! 10
    end if
  end do

  do i = 1, n_prob-1
    wall_p(i) = dens_const * wall_p(i)*u_tau_p/visc_const
    u_p   (i) = u_p(i) / u_tau_p
    v_p   (i) = v_p(i) / u_tau_p
    w_p   (i) = w_p(i) / u_tau_p

    uu_p (i) = uu_p (i) / (u_tau_p**2)
    vv_p (i) = vv_p (i) / (u_tau_p**2)
    ww_p (i) = ww_p (i) / (u_tau_p**2)
    uw_p (i) = uw_p (i) / (u_tau_p**2)

    kin_p   (i) = kin_p(i) / u_tau_p**2                           ! kin % n(c)
    eps_p   (i) = eps_p(i)*visc_const / (u_tau_p**4*dens_const)   ! eps % n(c)
    uw_mod_p(i) = uw_mod_p (i) / (u_tau_p**2)
    ww_mod_p(i) = ww_mod_p (i) / (u_tau_p**2)
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(fu2,'(10e15.7)')                                         &
                        wall_p(i),                                   &  !  1
                        u_p(i),                                      &  !  2
                        0.5*(uu_p(i)+vv_p(i)+ww_p(i)),               &  !  3
                        kin_p(i),                                    &  !  4
                        (0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i)),  &  !  5
                        uw_p(i),                                     &  !  6
                        uw_mod_p(i),                                 &  !  7
                        (uw_p(i) + uw_mod_p(i)),                     &  !  8
                        vis_t_p(i),                                  &  !  9
                        (ww_p(i) + ww_mod_p(i))                         ! 10
    end if
  end do

  close(fu1)
  close(fu2)

  if(First_Proc())  print *, '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
