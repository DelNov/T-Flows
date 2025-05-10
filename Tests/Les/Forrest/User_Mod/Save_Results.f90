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
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, optional        :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Face_Type), pointer :: flux
  integer                  :: fu, n_prob, pl, c, i, count, n_points, n_stat
  character(SL)            :: coord_name, res_name
  real, allocatable        :: z_p(:), ind(:), wall_p(:),                 &
                              u_p(:), v_p(:), w_p(:), y_plus_p(:),       &
                              kin_p(:), eps_p(:), uw_p(:), uw_mod_p(:),  &
                              uu_p(:), vv_p(:), ww_p(:), lai_p(:),     &
                              t_p(:), t2_p(:), ut_p(:), vt_p(:), wt_p(:)
  integer,allocatable      :: n_p(:), n_count(:)
  real                     :: visc_const
  real                     :: u_ref, uw_ref, h_ref
  logical                  :: there
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return
  if(.not. Turb % statistics) return

  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_TURB_STATISTICS',  &
                               HUGE_INT, n_stat, .false.)

  if(Time % Curr_Dt() < n_stat) return

  h_ref = 20.0

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  ! Take constant physical properties
  call Control % Dynamic_Viscosity(visc_const)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  ! Set file names for results
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
      print *, '# 3 0.2'
      print *, '# ... '
      print *, '#--------------------------------------------------------------'
    end if

    return
  end if

  n_points = 0

  call File % Open_For_Reading_Ascii(coord_name, fu)
  open(fu, file=coord_name)

  ! Write the number of searching intervals
  read(fu,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(fu,*) ind(pl), z_p(pl)
  end do
  close(fu)

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
  allocate(lai_p   (n_prob));  lai_p    = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0
  if(Flow % heat_transfer) then
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
    do c = Cells_In_Domain()
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

        lai_p (i) = lai_p (i)                           &
             +     3.56426355e-6 * min(Grid % wall_dist(c),20.0)**5       &
                 - 1.35550889e-4 * min(Grid % wall_dist(c),20.0)**4       &
                 + 1.71259681e-4 * min(Grid % wall_dist(c),20.0)**3       &
                 + 2.08486542e-2 * min(Grid % wall_dist(c),20.0)**2       &
                 + 4.22349842e-3 * min(Grid % wall_dist(c),20.0)          &
                 + 5.04175688e-1

        if(Turb % model .eq. HYBRID_LES_RANS) then 
          kin_p   (i) = kin_p   (i) + Turb % kin_mean(c)
          eps_p   (i) = eps_p   (i) + Turb % eps_mean(c)
          uw_mod_p(i) = uw_mod_p(i) + Turb % vis_t_eff(c)*(u % z(c) + w % x(c))
        end if

        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)

        if(Flow % heat_transfer) then
          t_p (i) = t_p (i) + Turb % t_mean(c)
          if(Turb % model .eq. HYBRID_LES_RANS) then 
            t2_p(i) = t2_p(i) + Turb % t2_mean(c)  &
                              - Turb % t_mean(c) * Turb % t_mean(c)
            ut_p(i) = ut_p(i) + Turb % ut_res(c)   &
                              - Turb % u_mean(c) * Turb % t_mean(c)
            vt_p(i) = vt_p(i) + Turb % vt_res(c)   &
                              - Turb % v_mean(c) * Turb % t_mean(c)
            wt_p(i) = wt_p(i) + Turb % wt_res(c)   &
                              - Turb % w_mean(c) * Turb % t_mean(c)
          end if
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

    call Global % Sum_Real(uu_p    (pl))
    call Global % Sum_Real(vv_p    (pl))
    call Global % Sum_Real(ww_p    (pl))
    call Global % Sum_Real(uw_p    (pl))
    call Global % Sum_Real(uw_mod_p(pl))
    call Global % Sum_Real(kin_p   (pl))
    call Global % Sum_Real(eps_p   (pl))
    call Global % Sum_Real(lai_p   (pl))
    call Global % Sum_Real(y_plus_p(pl))

    count =  count + n_count(pl)

    if(Flow % heat_transfer) then
      call Global % Sum_Real(t_p (pl))
      call Global % Sum_Real(t2_p(pl))
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
      uu_p  (i) = uu_p  (i) / n_count(i)
      vv_p  (i) = vv_p  (i) / n_count(i)
      ww_p  (i) = ww_p  (i) / n_count(i)
      uw_p  (i) = uw_p  (i) / n_count(i)

      uw_mod_p(i) = uw_mod_p(i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      lai_p (i) = lai_p (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
      if(Flow % heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        t2_p(i) = t2_p(i) / n_count(i)
        ut_p(i) = ut_p(i) / n_count(i)
        vt_p(i) = vt_p(i) / n_count(i)
        wt_p(i) = wt_p(i) / n_count(i)
      end if
    end if
  end do

  u_ref  = 0.0
  uw_ref = 0.0 
  do i = 1, n_prob
    u_ref  = max(u_ref,u_p(i))
    uw_ref = max(uw_ref,abs(uw_p(i)))
  end do

  call File % Open_For_Writing_Ascii(res_name, fu)
  open(fu, file=res_name)

  if(Turb % model .eq. HYBRID_LES_RANS) then
    if(Flow % heat_transfer) then
      write(fu,'(a1,2x,a105)') '#',' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved, kin_modeled,' //  &
                                   ' kin_tot,  uw_resolved,'     //  &
                                   ' uw_modeled, uw_tot, lai'    //  &
                                    ' t mean, ut_res, vt_res, wt_res,'
      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(14es15.5e3)') wall_p(i)/h_ref,                         &
                                  u_p(i)/u_ref,                             &
                                  0.5*(uu_p(i)+vv_p(i)+ww_p(i))/u_ref**2,   &
                                  kin_p(i)/u_ref**2,                        &
                                  (0.5*(uu_p(i)+vv_p(i)+ww_p(i))            &
                                  + kin_p(i))/u_ref**2 ,                    &
                                  abs(uw_p(i))/uw_ref,                      &
                                  abs(uw_mod_p(i))/uw_ref,                  &
                                  (abs(uw_p(i)) + abs(uw_mod_p(i)))/uw_ref, &
                                  lai_p(i),                                 &
                                  t_p(i),                                   &
                                  t2_p(i),                                  &
                                  ut_p(i),                                  &
                                  vt_p(i),                                  &
                                  wt_p(i)
        end if
      end do
    else
      write(fu,'(a1,2x,a85)') '#', ' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved, kin_modeled,' //  &
                                   ' kin_tot,  uw_resolved,'     //  &
                                   ' uw_modeled, uw_tot, lai'
      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(10es15.5e3)') wall_p(i)/h_ref,                         &
                                  u_p(i)/u_ref,                             &
                                  0.5*(uu_p(i)+vv_p(i)+ww_p(i))/u_ref**2,   &
                                  kin_p(i)/u_ref**2,                        &
                                  0.5*(uu_p(i)+vv_p(i)+ww_p(i)              &
                                  + kin_p(i))/u_ref**2,                     &
                                  abs(uw_p(i))/uw_ref,                      &
                                  abs(uw_mod_p(i))/uw_ref,                  &
                                  (abs(uw_p(i)) + abs(uw_mod_p(i)))/uw_ref, &
                                  lai_p(i) 
        end if
      end do
    end if
  else if(Turb % model .eq. LES_DYNAMIC) then
    if(Flow % heat_transfer) then
      write(fu,'(a1,2x,a105)') '#',' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved,'              //  &
                                   ' uw_resolved,'               //  &
                                   ' lai'                        //  &
                                   ' t mean, ut_res, vt_res, wt_res,'
      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(10es15.5e3)') wall_p(i)/h_ref,                         &
                                  u_p(i)/u_ref,                             &
                                  0.5*(uu_p(i)+vv_p(i)+ww_p(i))/u_ref**2,   &
                                  abs(uw_p(i))/uw_ref,                      &
                                  lai_p(i),                                 &
                                  t_p(i),                                   &
                                  t2_p(i),                                  &
                                  ut_p(i),                                  &
                                  vt_p(i),                                  &
                                  wt_p(i)
        end if
      end do
    else
      write(fu,'(a1,2x,a55)') '#', ' z,'                         //  &
                                   ' u mean,'                    //  &
                                   ' kin_resolved,'              //  &
                                   ' uw_resolved,'               //  &
                                   ' lai'
      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(5es15.5e3)') wall_p(i)/h_ref,                           &
                                  u_p(i)/u_ref,                              &
                                  0.5*(uu_p(i)+vv_p(i)+ww_p(i))/u_ref**2,    &
                                  abs(uw_p(i))/uw_ref,                       &
                                  lai_p(i)
        end if
      end do
    end if
  end if

  close(fu)

  if(First_Proc())  print '(a)', ' # Finished with User_Mod_Save_Results.f90.'

  end subroutine
