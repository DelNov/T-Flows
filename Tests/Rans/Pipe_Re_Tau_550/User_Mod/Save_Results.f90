!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1r file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
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
  type(Var_Type),  pointer :: u, v, w, t
  type(Var_Type),  pointer :: kin, eps, zeta, f22
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: ut, vt, wt
  integer                  :: n_prob, pl, c, i, count, s, c1, c2, n_points, reg
  character(SL)            :: coord_name, res_name, res_name_plus
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:), y_plus_p(:),  &
                              kin_p(:), eps_p(:), f22_p(:),                 &
                              vis_t_p(:), uw_p(:), zeta_p(:),               &
                              t_p(:), tt_p(:), ut_p(:), vt_p(:), wt_p(:),   &
                              ind(:), wall_p(:)
  real, pointer            :: vis_t(:)
  integer,allocatable      :: n_p(:), n_count(:)
  real                     :: t_wall, t_tau, d_wall, nu_mean, b11, b12, rad, r
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p, t_inf
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  logical                  :: there
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  vis_t  => Turb % vis_t
  call Flow % Alias_Momentum    (u, v, w)
  call Flow % Alias_Energy      (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)

  ! Take constant physical properties
  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1r')

  ! Set file name for results
  call File % Set_Name(res_name,                      &
                       time_step = Time % Curr_Dt(),  &
                       appendix  = '-res',            &
                       extension = '.dat')
  call File % Set_Name(res_name_plus,                 &
                       time_step = Time % Curr_Dt(),  &
                       appendix  = '-res-plus',       &
                       extension = '.dat')

  !------------------!
  !   Read 1r file   !
  !------------------!
  inquire(file=coord_name, exist=there)
  if(.not. there) then
    if(First_Proc()) then
      print *, '#=============================================================='
      print *, '# In order to extract profiles and write them in ascii files'
      print *, '# the code has to read cell-faces coordinates '
      print *, '# in wall-normal direction in the ascii file ''case_name.1r.'''
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

  ubulk    = bulk % flux_z / bulk % area_z
  t_wall   = 0.0
  nu_mean  = 0.0
  n_points = 0

  open(9, file=coord_name)

  ! Write the number of searching intervals
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl = 1, n_prob
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

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(Flow % heat_transfer) then
    allocate(t_p (n_prob));  t_p  = 0.0
    allocate(tt_p(n_prob));  tt_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = Cells_In_Domain()
      rad = 1.0 - Grid % wall_dist(c)
      if( rad < (z_p(i)) .and.  &
          rad > (z_p(i+1))) then
        r = sqrt(Grid % xc(c)**2 + Grid % yc(c)**2)
        b11 = Grid % xc(c)/r
        b12 = Grid % yc(c)/r

        wall_p(i) = wall_p(i) + Grid % wall_dist(c)
        u_p(i)   = u_p(i) + u % n(c)
        v_p(i)   = v_p(i) + v % n(c)
        w_p(i)   = w_p(i) + w % n(c)

        kin_p(i)   = kin_p(i) + kin % n(c)
        eps_p(i)   = eps_p(i) + eps % n(c)
        uw_p(i)    = uw_p(i)  &
                   + b11 * vis_t(c) *(u % z(c) + w % x(c)) &
                   + b12 * vis_t(c) *(v % z(c) + w % y(c))
        vis_t_p(i) = vis_t_p(i) + vis_t(c) / visc_const
        y_plus_p(i)= y_plus_p(i) + Turb % y_plus(c)

        if(Turb % model .eq. K_EPS_ZETA_F) then
          f22_p(i)  = f22_p(i) + f22  % n(c)
          zeta_p(i) = zeta_p(i) + zeta % n(c)
        end if

        if(Flow % heat_transfer) then
          t_p(i)    = t_p(i)  + t % n(c)
          ut_p(i)   = ut_p(i) + ut % n(c)
          vt_p(i)   = vt_p(i) + vt % n(c)
          wt_p(i)   = wt_p(i) + wt % n(c)
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

    call Global % Sum_Real(kin_p(pl))
    call Global % Sum_Real(eps_p(pl))
    call Global % Sum_Real(uw_p(pl))
    call Global % Sum_Real(vis_t_p(pl))
    call Global % Sum_Real(y_plus_p(pl))

    call Global % Sum_Real(f22_p(pl))
    call Global % Sum_Real(zeta_p(pl))

    count =  count + n_count(pl)

    if(Flow % heat_transfer) then
      call Global % Sum_Real(t_p(pl))
      call Global % Sum_Real(tt_p(pl))
      call Global % Sum_Real(ut_p(pl))
      call Global % Sum_Real(vt_p(pl))
      call Global % Sum_Real(wt_p(pl))
    end if
  end do


  call Global % Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count(i)
      u_p   (i) = u_p   (i) / n_count(i)
      v_p   (i) = v_p   (i) / n_count(i)
      w_p   (i) = w_p   (i) / n_count(i)

      kin_p   (i) = kin_p(i)    / n_count(i)
      eps_p   (i) = eps_p(i)    / n_count(i)
      uw_p    (i) = uw_p(i)     / n_count(i)
      vis_t_p (i) = vis_t_p(i)  / n_count(i)
      f22_p   (i) = f22_p(i)    / n_count(i)
      zeta_p  (i) = zeta_p(i)   / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)

      if(Flow % heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        tt_p(i) = tt_p(i) / n_count(i)
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
                       abs(bulk % p_drop_z))/dens_const)
  else  
    u_tau_p =  sqrt( (visc_const*sqrt(u_p(1)**2 +   &
                                        v_p(1)**2 +   &
                                        w_p(1)**2)    &
                                        / wall_p(1))  &
                                        / dens_const)
  end if

  if(u_tau_p .eq. 0.0) then
    if(First_Proc()) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if

    return
  end if

  if(Flow % heat_transfer) then
    d_wall = 0.0
    do c = Cells_In_Domain()
      if(Grid % wall_dist(c) > d_wall) then
        d_wall = Grid % wall_dist(c)
        t_inf  = t % n(c)
      end if
    end do

    call Global % Wait

    if(Flow % heat_flux > 0.0) then
      call Global % Min_Real(t_inf)
    else
      call Global % Max_Real(t_inf)
    end if

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL .or.  &
         Grid % region % type(reg) .eq. WALLFL) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          t_wall   = t_wall + t % n(c2)
          nu_mean  = nu_mean + t % q(c2) / (cond_const*(t % n(c2) - t_inf))
          n_points = n_points + 1
        end do  ! faces
      end if    ! boundary condition type
    end do      ! regions

    call Global % Sum_Real(t_wall)
    call Global % Sum_Real(nu_mean)
    call Global % Sum_Int(n_points)

    call Global % Wait

    t_wall  = t_wall / n_points
    nu_mean = nu_mean / n_points
    t_tau   = Flow % heat_flux / (dens_const * capa_const * u_tau_p)
  end if

  open(3, file = res_name)
  open(4, file = res_name_plus)

  do i = 3, 4
    pr = visc_const * capa_const / cond_const
    re = dens_const * ubulk * 2.0 / visc_const
    cf_dean = 0.0791*(re)**(-0.25)
    cf      = u_tau_p**2/(0.5*ubulk**2)
    error   = abs(cf_dean - cf)/cf_dean * 100.0
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Re       = ', dens_const * ubulk * 2.0/visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', dens_const*u_tau_p/visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(i,'(a1,(a12,F12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(i,'(a1,(a12,f12.6,a2,a22))') & 
    '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'
    if(Flow % heat_transfer) then
      write(i,'(a1,(a12, F12.6))')'#', 'Nu number =', nu_mean 
      write(i,'(a1,(a12, F12.6,a2,A39))')'#', 'Nu_error  =',  &
            abs(0.023*0.5*re**0.8*pr**0.4 - nu_mean)          & 
            / (0.023*0.5*re**0.8*pr**0.4) * 100.0, ' %',      &
      'correlation of Dittus-Boelter is used.' 
    end if

    if(Turb % model .eq. K_EPS) then
      if(Flow % heat_transfer) then
        write(i,'(a1,2x,a60)') '#',  ' r,'                    //  &    !  1
                                     ' w,'                    //  &    !  2
                                     ' kin, eps, uw,'         //  &    !  3, 4, 5
                                     ' vis_t/visc_const,'     //  &    !  6
                                     ' t, ut, vt, wt,'                 !  7 - 10
      else
        write(i,'(a1,2X,a60)') '#', ' r,'                     //  &    !  1
                                    ' w,'                     //  &    !  2
                                    ' kin, eps, uw, vis_t/visc_const'  !  3-6
      end if
    else if(Turb % model .eq. K_EPS_ZETA_F) then
      if(Flow % heat_transfer) then
        write(i,'(a1,2x,a60)') '#',  ' r,'                    //  &  !  1
                                     ' w,'                    //  &  !  2
                                     ' kin, eps, uw,'         //  &  !  3, 4, 5
                                     ' f22, zeta,'            //  &  !  6, 7
                                     ' vis_t/visc_const,'     //  &  !  8 - 11
                                     ' t, ut, vt, wt'
      else
        write(i,'(a1,2x,a50)') '#', ' r,'                     //  &  !  1
                                    ' w,'                     //  &  !  2
                                    ' kin, eps, uw,'          //  &  !  3, 4, 5
                                    ' f22, zeta'              //  &  !  6, 7
                                    ' vis_t/visc_const,'             !  8
      end if
    end if
  end do

  if(Flow % heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(12es15.5e3,i5)') wall_p(i),   &  !  1
                                   w_p(i),      &  !  2
                                   kin_p(i),    &  !  3
                                   eps_p(i),    &  !  4
                                   uw_p(i),     &  !  5
                                   f22_p(i),    &  !  6
                                   zeta_p(i),   &  !  7
                                   vis_t_p(i),  &  !  8
                                   t_p(i),      &  !  9
                                   ut_p(i),     &  ! 10
                                   vt_p(i),     &  ! 11
                                   wt_p(i),     &  ! 12
                                   n_count(i)      ! 13
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(8es15.5e3,i5)')  wall_p(i),   &  !  1
                                   w_p(i),      &  !  2
                                   kin_p(i),    &  !  3
                                   eps_p(i),    &  !  4
                                   uw_p(i),     &  !  5
                                   f22_p(i),    &  !  6
                                   zeta_p(i),   &  !  7
                                   vis_t_p(i),  &  !  8
                                   n_count(i)      !  9
      end if
    end do
  end if

  close(3)

  do i = 1, n_prob-1
    wall_p(i) = dens_const * wall_p(i)*u_tau_p/visc_const
    w_p   (i) = w_p  (i) / u_tau_p
    kin_p (i) = kin_p(i) / u_tau_p**2                      ! kin%n(c)
    eps_p (i) = eps_p(i)*visc_const / (u_tau_p**4*dens_const)  ! eps%n(c)
    uw_p  (i) = uw_p (i) / (u_tau_p**2*dens_const)    ! vis_t(c)*(u%z(c)+w%x(c))

    if(Turb % model .eq. K_EPS_ZETA_F) then
      f22_p(i) = f22_p(i) * visc_const / u_tau_p**2  ! f22%n(c)
    end if

    if(Flow % heat_transfer) then
      t_p (i) = (t_wall - t_p(i)) / t_tau   ! t % n(c)
      ut_p(i) = ut_p(i) / (u_tau_p*t_tau)   ! ut % n(c)
      vt_p(i) = vt_p(i) / (u_tau_p*t_tau)   ! vt % n(c)
      wt_p(i) = wt_p(i) / (u_tau_p*t_tau)   ! wt % n(c)
    end if
  end do

  if(Flow % heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(12es15.5e3)') wall_p(i),   &  !  1
                                w_p(i),      &  !  2
                                kin_p(i),    &  !  3
                                eps_p(i),    &  !  4
                                uw_p(i),     &  !  5
                                f22_p(i),    &  !  6
                                zeta_p(i),   &  !  7
                                vis_t_p(i),  &  !  8
                                t_p(i),      &  !  9
                                ut_p(i),     &  ! 10
                                vt_p(i),     &  ! 11
                                wt_p(i)         ! 12
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(12es15.5e3)') wall_p(i),   &  !  1
                                w_p(i),      &  !  2
                                kin_p(i),    &  !  3
                                eps_p(i),    &  !  4
                                uw_p(i),     &  !  5
                                f22_p(i),    &  !  6
                                zeta_p(i),   &  !  7
                                vis_t_p(i)      !  8
      end if
    end do
  end if

  close(4)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(kin_p)
  deallocate(eps_p)
  deallocate(uw_p)
  deallocate(vis_t_p)
  deallocate(f22_p)
  deallocate(zeta_p)
  deallocate(y_plus_p)
  if(Flow % heat_transfer) then
    deallocate(t_p)
    deallocate(tt_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(First_Proc())  print *, '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
