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
  type(Var_Type),  pointer :: kin, eps, zeta, f, ut, vt, wt, t2
  integer             :: n_prob, pl, c, i, count, s, c1, c2, n_points, fu
  character(len=SL)   :: coord_name, res_name, res_name_plus, nu_name
  real, allocatable   :: z_p(:), tz_p(:), ti_p(:), w_p(:), t_p(:),         &
                         y_plus_p(:),  ind(:),  wall_p(:), kin_p(:),       &
                         eps_p(:), kin_mod_p(:),                           &
                         uw_p(:), uu_p(:), vv_p(:), ww_p(:),               &
                         t2_p(:), ut_p(:), vt_p(:), wt_p(:), t2_mod_p(:),  &
                         ut_mod(:), vt_mod(:), wt_mod(:), var_1(:),        &
                         var_2(:), var_3(:)
  integer, allocatable :: n_p(:), n_count(:)
  real                 :: t_wall, t_tau, d_wall, t_hot, t_cold, t_diff, pr_t
  real                 :: vel_ref
  logical              :: there
  real, contiguous, pointer :: tz_mean(:)
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  call Work % Connect_Real_Cell(tz_mean)

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)
  call Turb % Alias_T2          (t2)

  ! Set some constants
  t_cold =   5.0
  t_hot  =  15.0
  t_diff =  t_hot - t_cold
!  vel_ref= sqrt(t_diff*abs(grav_z)*Flow % beta)
  vel_ref = 21.28e-6

  call Flow % Grad_Component(Grid, Turb % t_mean, 3, tz_mean)
  call Flow % Grad_Variable(t)
  call Flow % Grad_Variable(Flow % w)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  ! Set file names for results
  call File % Set_Name(res_name,                      &
                       time_step = Time % Curr_Dt(),  &
                       appendix  = '-res',            &
                       extension = '.dat')
  call File % Set_Name(res_name_plus,                 &
                       time_step = Time % Curr_Dt(),  &
                       appendix  = '-res-plus',       &
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

  t_wall   = 0.0
  n_points = 0

  ! This functions opens a file in first available unit,
  ! and stores it in variable "fu".  Safer like this,
  ! you can't hit a file unit already open
  call File % Open_For_Reading_Ascii(res_name, fu)

  ! Write the number of searching intervals
  read(fu,*)  n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(fu, *) ind(pl), z_p(pl)
  end do
  close(fu)

  allocate(n_p      (n_prob));  n_p       = 0
  allocate(wall_p   (n_prob));  wall_p    = 0.0
  allocate(tz_p     (n_prob));  tz_p      = 0.0  ! dT / dz
  allocate(ti_p     (n_prob));  ti_p      = 0.0  ! T instant
  allocate(w_p      (n_prob));  w_p       = 0.0
  allocate(uu_p     (n_prob));  uu_p      = 0.0
  allocate(vv_p     (n_prob));  vv_p      = 0.0
  allocate(ww_p     (n_prob));  ww_p      = 0.0
  allocate(uw_p     (n_prob));  uw_p      = 0.0
  allocate(kin_mod_p(n_prob));  kin_mod_p = 0.0
  allocate(kin_p    (n_prob));  kin_p     = 0.0
  allocate(var_1    (n_prob));  var_1     = 0.0
  allocate(var_2    (n_prob));  var_2     = 0.0
  allocate(var_3    (n_prob));  var_3     = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(Flow % heat_transfer) then
    allocate(t_p (n_prob));     t_p      = 0.0
    allocate(t2_p(n_prob));     t2_p     = 0.0
    allocate(t2_mod_p(n_prob)); t2_mod_p = 0.0
    allocate(ut_p(n_prob));     ut_p     = 0.0
    allocate(vt_p(n_prob));     vt_p     = 0.0
    allocate(wt_p(n_prob));     wt_p     = 0.0
    allocate(ut_mod(n_prob));   ut_mod   = 0.0
    allocate(vt_mod(n_prob));   vt_mod   = 0.0
    allocate(wt_mod(n_prob));   wt_mod   = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
      pr_t = Turb % Prandtl_Turb(c)
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + Grid % zc(c)
        tz_p  (i) = tz_p  (i) - Flow % conductivity(c) * tz_mean(c)
        ti_p  (i) = ti_p  (i) + Turb % t_mean(c)
        w_p   (i) = w_p   (i) + Turb % w_mean(c)

        uu_p(i)   = uu_p(i) + Turb % uu_res(c)  &
                            - Turb % u_mean(c) * Turb % u_mean(c)
        vv_p(i)   = vv_p(i) + Turb % vv_res(c)  &
                            - Turb % v_mean(c) * Turb % v_mean(c)
        ww_p(i)   = ww_p(i) + Turb % ww_res(c)  &
                            - Turb % w_mean(c) * Turb % w_mean(c)
        uw_p(i)   = uw_p(i) + Turb % uw_res(c)  &
                            - Turb % u_mean(c) * Turb % w_mean(c)
        kin_p(i)  = kin_p(i) &
                + 0.5*(Turb % uu_res(c) - Turb % u_mean(c) * Turb % u_mean(c) &
                     + Turb % vv_res(c) - Turb % v_mean(c) * Turb % v_mean(c) &
                     + Turb % ww_res(c) - Turb % w_mean(c) * Turb % w_mean(c))
        kin_mod_p(i) = kin_mod_p(i) + Turb % kin_mean(c)

        if(Flow % heat_transfer) then
          t_p(i)  = t_p(i)  + (Turb % t_mean(c) - t_cold)/t_diff
          ut_p(i) = ut_p(i) + Turb % ut_res(c)  &
                            - Turb % u_mean(c) * Turb % t_mean(c)
          vt_p(i) = vt_p(i) + Turb % vt_res(c)  &
                            - Turb % v_mean(c) * Turb % t_mean(c)
          wt_p(i) = wt_p(i) + Turb % wt_res(c)  &
                            - Turb % w_mean(c) * Turb % t_mean(c)
          t2_p(i) = t2_p(i) + Turb % t2_res(c)  &
                            - Turb % t_mean(c) * Turb % t_mean(c)
!         var_3(i) = var_3(i) - Flow % conductivity(c)* t % z(c) &
!                             + Turb % wt_res(c) + Turb % wt % n(c)
!         if(Turb % model .eq. HYBRID_LES_RANS) then
            ut_mod(i)   = ut_mod(i)   + Turb % ut % n(c)
            vt_mod(i)   = vt_mod(i)   + Turb % vt % n(c)
            wt_mod(i)   = wt_mod(i)   + Turb % wt_mean(c)
            t2_mod_p(i) = t2_mod_p(i) + Turb % t2_mean(c)
            var_1(i)    = var_1(i)   &
                        + 0.2 * Turb % t_scale(c)              &
                       ! * (0.666*Turb % kin_mean(c)*tz_mean(c)  &  ! t % z(c)
                              * (tz_mean(c) &                       ! t % z(c)
                        + 0.0*0.6*Turb % wt % n(c) * w % z(c)  &
                        + 0.0*0.6*0.0327*Turb % t2_mean(c))
            var_2(i)  = var_2(i) + 0.2 * Turb % t_scale(c)  &
                                 * 0.6*0.0327 * Turb % t2_mean(c)
                               !   0.2 * Turb % t_scale(c)   &
                               ! * 0.666*0.6*0.0327*Turb % t2_mean(c)
                               !   Turb % kin % n(c) * t % z(c)
!         end if
          var_3(i) = var_3(i)                 &
                   + 0.2 * Turb % t_scale(c)  &
                         * (0.666*Turb % kin_mean(c)*tz_mean(c))    ! t % z(c)
        end if
        n_count(i) = n_count(i) + 1
      end if
    end do
  end do


  ! Average over all processors
  do pl=1, n_prob-1
    call Global % Sum_Int(n_count(pl))

    call Global % Sum_Real(wall_p(pl))

    call Global % Sum_Real(tz_p(pl))
    call Global % Sum_Real(ti_p(pl))
    call Global % Sum_Real(w_p(pl))

    call Global % Sum_Real(uu_p(pl))
    call Global % Sum_Real(vv_p(pl))
    call Global % Sum_Real(ww_p(pl))
    call Global % Sum_Real(uw_p(pl))
    call Global % Sum_Real(kin_p(pl))
    call Global % Sum_Real(kin_mod_p(pl))
    call Global % Sum_Real(var_1(pl))
    call Global % Sum_Real(var_2(pl))
    call Global % Sum_Real(var_3(pl))

    count =  count + n_count(pl)

    if(Flow % heat_transfer) then
      call Global % Sum_Real(t_p(pl))
      call Global % Sum_Real(t2_p(pl))
      call Global % Sum_Real(t2_mod_p(pl))
      call Global % Sum_Real(ut_p(pl))
      call Global % Sum_Real(vt_p(pl))
      call Global % Sum_Real(wt_p(pl))
      call Global % Sum_Real(ut_mod(pl))
      call Global % Sum_Real(vt_mod(pl))
      call Global % Sum_Real(wt_mod(pl))
    end if
  end do

  call Global % Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count(i)
      tz_p  (i) = tz_p (i)  / n_count(i)
      ti_p  (i) = ti_p (i)  / n_count(i)
      w_p   (i) = w_p  (i)  / n_count(i)

      uu_p(i)      = uu_p(i)      / n_count(i)
      vv_p(i)      = vv_p(i)      / n_count(i)
      ww_p(i)      = ww_p(i)      / n_count(i)
      uw_p(i)      = uw_p(i)      / n_count(i)
      kin_p(i)     = kin_p(i)     / n_count(i)
      kin_mod_p(i) = kin_mod_p(i) / n_count(i)
      var_1(i)     = var_1(i)     / n_count(i)
      var_2(i)     = var_2(i)     / n_count(i)
      var_3(i)     = var_3(i)     / n_count(i)

      if(Flow % heat_transfer) then
        t_p (i)     = t_p (i)     / n_count(i)
        t2_p(i)     = t2_p(i)     / n_count(i)
        t2_mod_p(i) = t2_mod_p(i) / n_count(i)
        ut_p(i)     = ut_p(i)     / n_count(i)
        vt_p(i)     = vt_p(i)     / n_count(i)
        wt_p(i)     = wt_p(i)     / n_count(i)
        ut_mod(i)   = ut_mod(i)   / n_count(i)
        vt_mod(i)   = vt_mod(i)   / n_count(i)
        wt_mod(i)   = wt_mod(i)   / n_count(i)
      end if
    end if
  end do

  ! This functions opens a file in first available unit,
  ! and stores it in variable "fu".  Safer like this,
  ! you can't hit a file unit already open
  call File % Open_For_Writing_Ascii(res_name, fu)

  if(Flow % heat_transfer) then
    if(First_Proc()) then
      write(*,'(a1,(a12, f12.6))')'#', ' Nu number = ',  &
               tz_p(1) / (t_hot - t_cold)
    end if
  end if

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(fu,'(15es15.5e3)')                                 &
            wall_p(i),                                         &  !  1
            (ti_p(i) - t_cold)/t_diff,                         &  !  2
            sqrt(kin_p(i))/vel_ref,                            &  !  3
            sqrt(kin_mod_p(i))/vel_ref,                        &  !  4
            sqrt(kin_p(i) + kin_mod_p(i))/vel_ref,             &  !  5
            sqrt(t2_p(i))/10.0,                                &  !  6
            sqrt(t2_mod_p(i))/10.0,                            &  !  7
            sqrt(t2_p(i)+t2_mod_p(i))/10.0,                    &  !  8
            wt_p(i)/(vel_ref*10.0),                            &  !  9
            wt_mod(i)/(vel_ref*10.0),                          &  ! 10
            tz_p(i)/(vel_ref*10.0),                            &  ! 11
            (wt_p(i) + wt_mod(i) + tz_p(i))/(vel_ref*10.0),    &  ! 12
            var_1(i)/(vel_ref*10.0), var_2(i)/(vel_ref*10.0),  &  ! 13
            var_3(i)/(vel_ref*10.0)                               ! 14
    end if
  end do

  close(fu)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(tz_p)
  deallocate(ti_p)
  deallocate(w_p)
  deallocate(uu_p)
  deallocate(vv_p)
  deallocate(ww_p)
  deallocate(uw_p)
  deallocate(kin_p)
  deallocate(kin_mod_p)
  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  if(Flow % heat_transfer) then
    deallocate(t_p)
    deallocate(t2_p)
    deallocate(t2_mod_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(First_Proc())  print *, '# Finished with User_Mod_Save_Results.f90.'

  call Work % Disconnect_Real_Cell(tz_mean)

  end subroutine
