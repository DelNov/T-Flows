!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, swarm, ts)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: ts   ! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Face_Type), pointer :: flux
  integer                  :: n_prob, pl, c, i, count, s, c1, c2, n_points
  character(len=80)        :: coord_name, res_name, res_name_plus
  real, allocatable        :: z_p(:), ind(:), wall_p(:),                 &
                              u_p(:), v_p(:), w_p(:), y_plus_p(:),       &
                              kin_p(:), eps_p(:), uw_p(:), uw_mod_p(:),  &
                              uu_p(:), vv_p(:), ww_p(:), vis_t_p(:),     &
                              t_p(:), t2_p(:), ut_p(:), vt_p(:), wt_p(:)
  integer,allocatable      :: n_p(:), n_count(:)
  real                     :: t_wall, t_tau, d_wall, nu_mean, t_inf
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  logical                  :: there
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Field_Mod_Alias_Energy  (flow, t)

  ! Take constant physical properties
  call Control_Mod_Mass_Density        (dens_const)
  call Control_Mod_Dynamic_Viscosity   (visc_const)
  call Control_Mod_Heat_Capacity       (capa_const)
  call Control_Mod_Thermal_Conductivity(cond_const)

  ! Set the name for coordinate file
  call File_Mod_Set_Name(coord_name, extension='.1d')

  ! Set file names for results
  call File_Mod_Set_Name(res_name,         &
                         time_step=ts,     &   
                         appendix='-res',  &
                         extension='.dat')
  call File_Mod_Set_Name(res_name_plus,         &
                         time_step=ts,          &
                         appendix='-res-plus',  &
                         extension='.dat')

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file=coord_name, exist=there)
  if(.not. there) then
    if(this_proc < 2) then
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

  ubulk    = bulk % flux_x / (dens_const * bulk % area_x)
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
  allocate(uu_p    (n_prob));  uu_p     = 0.0
  allocate(vv_p    (n_prob));  vv_p     = 0.0
  allocate(ww_p    (n_prob));  ww_p     = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(uw_mod_p(n_prob));  uw_mod_p = 0.0
  allocate(vis_t_p (n_prob));  vis_t_p  = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0
  if(heat_transfer) then
    allocate(t_p (n_prob));  t_p  = 0.0
    allocate(t2_p(n_prob));  t2_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells 
      if(grid % zc(c) > (z_p(i)) .and.  &
         grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + grid % wall_dist(c)
        u_p   (i) = u_p   (i) + turb % u_mean(c)
        v_p   (i) = v_p   (i) + turb % v_mean(c)
        w_p   (i) = w_p   (i) + turb % w_mean(c)

        uu_p(i) = uu_p(i) + turb % uu_res(c)  &
                          - turb % u_mean(c) * turb % u_mean(c)
        vv_p(i) = vv_p(i) + turb % vv_res(c)  &
                          - turb % v_mean(c) * turb % v_mean(c)
        ww_p(i) = ww_p(i) + turb % ww_res(c)  &
                          - turb % w_mean(c) * turb % w_mean(c)
        uw_p(i) = uw_p(i) + turb % uw_res(c)  &
                          - turb % u_mean(c) * turb % w_mean(c)

        kin_p   (i) = kin_p   (i) + turb % kin_mean(c)
        eps_p   (i) = eps_p   (i) + turb % eps_mean(c)
        uw_mod_p(i) = uw_mod_p(i) + turb % vis_t_eff(c)*(u % z(c) + w % x(c))
        vis_t_p (i) = vis_t_p (i) + turb % vis_t_eff(c) / visc_const
        y_plus_p(i) = y_plus_p(i) + turb % y_plus(c)

        if(heat_transfer) then
          t_p (i) = t_p (i) + turb % t_mean(c)
          t2_p(i) = t2_p(i) + turb % t2_mean(c)  &
                            - turb % t_mean(c) * turb % t_mean(c)
          ut_p(i) = ut_p(i) + turb % ut_res(c)   &
                            - turb % u_mean(c) * turb % t_mean(c)
          vt_p(i) = vt_p(i) + turb % vt_res(c)   &
                            - turb % v_mean(c) * turb % t_mean(c)
          wt_p(i) = wt_p(i) + turb % wt_res(c)   &
                            - turb % w_mean(c) * turb % t_mean(c)
        end if
        n_count(i) = n_count(i) + 1
      end if
    end do
  end do

  ! Average over all processors
  do pl=1, n_prob-1
    call Comm_Mod_Global_Sum_Int(n_count(pl))

    call Comm_Mod_Global_Sum_Real(wall_p(pl))

    call Comm_Mod_Global_Sum_Real(u_p(pl))
    call Comm_Mod_Global_Sum_Real(v_p(pl))
    call Comm_Mod_Global_Sum_Real(w_p(pl))

    call Comm_Mod_Global_Sum_Real(uu_p    (pl))
    call Comm_Mod_Global_Sum_Real(vv_p    (pl))
    call Comm_Mod_Global_Sum_Real(ww_p    (pl))
    call Comm_Mod_Global_Sum_Real(uw_p    (pl))
    call Comm_Mod_Global_Sum_Real(uw_mod_p(pl))
    call Comm_Mod_Global_Sum_Real(kin_p   (pl))
    call Comm_Mod_Global_Sum_Real(eps_p   (pl))
    call Comm_Mod_Global_Sum_Real(vis_t_p (pl))
    call Comm_Mod_Global_Sum_Real(y_plus_p(pl))

    count =  count + n_count(pl)

    if(heat_transfer) then
      call Comm_Mod_Global_Sum_Real(t_p (pl))
      call Comm_Mod_Global_Sum_Real(t2_p(pl))
      call Comm_Mod_Global_Sum_Real(ut_p(pl))
      call Comm_Mod_Global_Sum_Real(vt_p(pl))
      call Comm_Mod_Global_Sum_Real(wt_p(pl))
    end if
  end do

  call Comm_Mod_Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count(i)
      u_p   (i) = u_p   (i) / n_count(i)
      v_p   (i) = v_p   (i) / n_count(i)
      w_p   (i) = w_p   (i) / n_count(i)
      uu_p  (i) = uu_p  (i) / n_count(i)
      vv_p  (i) = vv_p  (i) / n_count(i)
      ww_p  (i) = ww_p  (i) / n_count(i)
      uw_p  (i) = uw_p  (i) / n_count(i)

      uw_mod_p(i) = uw_mod_p(i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      vis_t_p (i) = vis_t_p (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
      if(heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        t2_p(i) = t2_p(i) / n_count(i)
        ut_p(i) = ut_p(i) / n_count(i)
        vt_p(i) = vt_p(i) / n_count(i)
        wt_p(i) = wt_p(i) / n_count(i)
      end if
    end if
  end do

  ! Calculating friction velocity and friction temperature
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
    if(this_proc < 2) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if
    return
  end if

  if(heat_transfer) then
    d_wall = 0.0 
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) > d_wall) then
        d_wall = grid % wall_dist(c)
        t_inf  = turb % t_mean(c)
      end if
    end do

    call Comm_Mod_Wait

    if(flow % heat_flux > 0.0) then
      call Comm_Mod_Global_Min_Real(t_inf)
    else
      call Comm_Mod_Global_Max_Real(t_inf)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if( Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. WALL .or.  &
            Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. WALLFL) then

          t_wall   = t_wall + turb % t_mean(c2)
          nu_mean  = nu_mean + t % q(c2)  &
                   / (cond_const*(turb % t_mean(c2) - t_inf))
          n_points = n_points + 1
        end if
      end if
    end do

    call Comm_Mod_Global_Sum_Real(t_wall)
    call Comm_Mod_Global_Sum_Real(nu_mean)
    call Comm_Mod_Global_Sum_Int(n_points)

    call Comm_Mod_Wait

    t_wall  = t_wall / n_points
    nu_mean = nu_mean / n_points
    t_tau   = flow % heat_flux / (dens_const * capa_const * u_tau_p)
  end if

  open(3, file = res_name)
  open(4, file = res_name_plus)

  do i = 3, 4
    pr = visc_const * capa_const / cond_const
    re = dens_const * ubulk * 2.0 / visc_const
    cf_dean = 0.073*(re)**(-0.25)
    cf      = u_tau_p**2/(0.5*ubulk**2)
    error   = abs(cf_dean - cf)/cf_dean * 100.0
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(i,'(a1,(a12,e12.6))')  &
    '#', 're       = ', dens_const * ubulk * 2.0 / visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', dens_const * u_tau_p / visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(i,'(a1,(a12,f12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(i,'(a1,(a12,f12.6,a2,a22))') & 
    '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'
    if(heat_transfer) then
      write(i,'(a1,(a12, f12.6))')'#', 'Nu number =', nu_mean
      write(i,'(a1,(a12, f12.6,a2,a39))')'#', 'Nu_error  =',  &
            abs(0.023*0.5*re**0.8*pr**0.4 - nu_mean)          &
            / (0.023*0.5*re**0.8*pr**0.4) * 100.0, ' %',      &
            'correlation of Dittus-Boelter is used.'
    end if

    if(heat_transfer) then
      write(i,'(a1,2X,a105)') '#', ' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved, kin_modeled,' //  &
                                   ' kin_tot,  uw_resolved,'     //  &
                                   ' uw_modeled, uw_tot, vis_t'  //  &
                                   ' t mean, ut_res, vt_res, wt_res,'
    else
      write(i,'(a1,2X,a85)') '#',  ' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved, kin_modeled,' //  &
                                   ' kin_tot,  uw_resolved,'     //  &
                                   ' uw_modeled, uw_tot, vis_t'
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(14es15.5e3)') wall_p(i),                                   &
                            u_p(i),                                          &
                            0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                   &
                            kin_p(i),                                        &
                            0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i),        &
                            abs(uw_p(i)),                                    &
                            abs(uw_mod_p(i)),                                &
                            abs(uw_p(i)) + abs(uw_mod_p(i)),                 &
                            vis_t_p(i),                                      &
                            t_p(i),                                          &
                            t2_p(i),                                         &
                            ut_p(i),                                         &
                            vt_p(i),                                         &
                            wt_p(i)
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(9es15.5e3)')  wall_p(i),                                   &
                                u_p(i),                                      &
                                0.5*(uu_p(i)+vv_p(i)+ww_p(i)),               &
                                kin_p(i),                                    &
                                (0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i)),  &
                                uw_p(i),                                     &
                                uw_mod_p(i),                                 &
                                (uw_p(i) + uw_mod_p(i)),                     &
                                vis_t_p(i)
      end if
    end do
  end if

  do i = 1, n_prob
    wall_p(i) = dens_const * wall_p(i)*u_tau_p / visc_const
    u_p   (i) = u_p(i) / u_tau_p
    v_p   (i) = v_p(i) / u_tau_p
    w_p   (i) = w_p(i) / u_tau_p

    kin_p (i) = kin_p(i) / u_tau_p**2                              ! kin % n(c)
    eps_p (i) = eps_p(i) * visc_const / (u_tau_p**4 * dens_const)  ! eps % n(c)
    uu_p  (i) = uu_p (i) / (u_tau_p**2)
    vv_p  (i) = vv_p (i) / (u_tau_p**2)
    ww_p  (i) = ww_p (i) / (u_tau_p**2)
    uw_p  (i) = uw_p (i) / (u_tau_p**2)
    uw_mod_p(i) = uw_mod_p (i) / (u_tau_p**2)

    if(heat_transfer) then
      t_p (i) = (t_wall - t_p(i)) / t_tau  ! t % n(c)
      t2_p(i) = t2_p(i) / (t_tau*t_tau)    ! ut % n(c)
      ut_p(i) = ut_p(i) / (u_tau_p*t_tau)  ! ut % n(c)
      vt_p(i) = vt_p(i) / (u_tau_p*t_tau)  ! vt % n(c)
      wt_p(i) = wt_p(i) / (u_tau_p*t_tau)  ! wt % n(c)
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(14es15.5e3)')                                & ! write
              wall_p(i),                                       & ! 1    
              u_p(i),                                          & ! 2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                   & ! 3
              kin_p(i),                                        & ! 4
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i),        & ! 5
              abs(uw_p(i)),                                    & ! 6
              abs(uw_mod_p(i)),                                & ! 6
              abs(uw_p(i)) + abs(uw_mod_p(i)),                 & ! 8
              vis_t_p(i),                                      & ! 9
              t_p(i),                                          & ! 10
              t2_p(i),                                         & ! 11 
              ut_p(i),                                         & ! 12
              vt_p(i),                                         & ! 13
              wt_p(i)                                            ! 14
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(9es15.5e3)')                             &  !  write
              wall_p(i),                                   &  !  1
              u_p(i),                                      &  !  2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),               &  !  3
              kin_p(i),                                    &  !  4
              (0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i)),  &  !  5
              uw_p(i),                                     &  !  6
              uw_mod_p(i),                                 &  !  7
              (uw_p(i) + uw_mod_p(i)),                     &  !  8
              vis_t_p(i)                                      !  9
      end if
    end do
  end if

  close(3)
  close(4)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(uu_p)
  deallocate(vv_p)
  deallocate(ww_p)
  deallocate(uw_p)
  deallocate(uw_mod_p)
  deallocate(kin_p)
  deallocate(eps_p)
  deallocate(vis_t_p)
  deallocate(y_plus_p)
  if(heat_transfer) then
    deallocate(t_p)
    deallocate(t2_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
