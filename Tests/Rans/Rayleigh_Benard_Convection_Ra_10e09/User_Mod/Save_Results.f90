!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, turb, Vof, swarm, ts)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: ts
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n_prob, pl, c, i, count, s, c1, c2, n_points
  character(SL)       :: coord_name, res_name, res_name_plus
  real, allocatable   :: z_p(:), tz_p(:), ti_p(:), w_p(:), t_p(:),         &
                         y_plus_p(:),  ind(:),  wall_p(:), kin_p(:),       &
                         eps_p(:), kin_mod_p(:),                           &
                         uw_p(:), uu_p(:), vv_p(:), ww_p(:),               &
                         t2_p(:), ut_p(:), vt_p(:), wt_p(:), t2_mod_p(:),  &
                         ut_mod(:), vt_mod(:), wt_mod(:)
  integer, allocatable :: n_p(:), n_count(:)
  real                 :: t_wall, t_tau, d_wall, t_hot, t_cold, t_diff
  logical              :: there
!==============================================================================!

!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: kin, eps, zeta, f, ut, vt, wt, t2
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f)
  call Turb_Mod_Alias_Stresses   (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt) 
  call Turb_Mod_Alias_T2          (turb, t2) 

  ! Set some constants
  t_cold =  5.0
  t_hot  = 20.0
  t_diff =  t_hot - t_cold

  call Flow % Grad_Variable(t)
  call Flow % Grad_Variable(Flow % u)
  call Flow % Grad_Variable(Flow % v)
  call Flow % Grad_Variable(Flow % w)

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  ! Set file names for results
  call File % Set_Name(res_name,            &
                       time_step = ts,      &
                       appendix  = '-res',  &
                       extension = '.dat')
  call File % Set_Name(res_name_plus,            &
                       time_step = ts,           &
                       appendix  = '-res-plus',  &
                       extension = '.dat')

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

  t_wall   = 0.0
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

  allocate(n_p      (n_prob));  n_p      = 0
  allocate(wall_p   (n_prob));  wall_p   = 0.0
  allocate(tz_p     (n_prob));  tz_p     = 0.0  ! dT / dz
  allocate(ti_p     (n_prob));  ti_p     = 0.0  ! T instant
  allocate(w_p      (n_prob));  w_p      = 0.0
  allocate(uu_p     (n_prob));  uu_p     = 0.0
  allocate(vv_p     (n_prob));  vv_p     = 0.0
  allocate(ww_p     (n_prob));  ww_p     = 0.0
  allocate(uw_p     (n_prob));  uw_p     = 0.0
  allocate(kin_mod_p(n_prob));  kin_mod_p= 0.0
  allocate(kin_p    (n_prob));  kin_p    = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(Flow % heat_transfer) then
    allocate(t_p (n_prob));     t_p = 0.0
    allocate(t2_p(n_prob));     t2_p = 0.0
    allocate(t2_mod_p(n_prob)); t2_mod_p = 0.0
    allocate(ut_p(n_prob));     ut_p = 0.0
    allocate(vt_p(n_prob));     vt_p = 0.0
    allocate(wt_p(n_prob));     wt_p = 0.0
    allocate(ut_mod(n_prob));   ut_mod = 0.0
    allocate(vt_mod(n_prob));   vt_mod = 0.0
    allocate(wt_mod(n_prob));   wt_mod = 0.0
    allocate(var_1(n_prob));    var_1  = 0.0
    allocate(var_2(n_prob));    var_2  = 0.0
    allocate(var_3(n_prob));    var_3  = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, Grid % n_cells - Grid % comm % n_buff_cells 
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + Grid % zc(c)
        tz_p  (i) = tz_p  (i) + t % z(c)
        ti_p  (i) = ti_p  (i) + t % n(c)
        w_p   (i) = w_p   (i) + turb % w_mean(c)

        uu_p(i)   = uu_p(i) + turb % uu_res(c)  &
                            - turb % u_mean(c) * turb % u_mean(c)
        vv_p(i)   = vv_p(i) + turb % vv_res(c)  &
                            - turb % v_mean(c) * turb % v_mean(c)
        ww_p(i)   = ww_p(i) + turb % ww_res(c)  &
                            - turb % w_mean(c) * turb % w_mean(c)
        uw_p(i)   = uw_p(i) + turb % uw_res(c)  &
                            - turb % u_mean(c) * turb % w_mean(c)
        kin_p(i)  = kin_p(i) &
                + 0.5*(turb % uu_res(c) - turb % u_mean(c) * turb % u_mean(c) &
                     + turb % vv_res(c) - turb % v_mean(c) * turb % v_mean(c) &
                     + turb % ww_res(c) - turb % w_mean(c) * turb % w_mean(c))
        kin_mod_p(i) = kin_mod_p(i) + turb % kin_mean(c)

        if(Flow % heat_transfer) then
          t_p(i)  = t_p(i)  + turb % t_mean(c) 
          t2_p(i) = t2_p(i) + turb % t2_res(c)  &
                            - turb % t_mean(c) * turb % t_mean(c)
          t2_mod_p(i) = t2_mod_p(i) + turb % t2_mean(c)
          ut_p(i) = ut_p(i) + turb % ut_res(c)  &
                            - turb % u_mean(c) * turb % t_mean(c)
          vt_p(i) = vt_p(i) + turb % vt_res(c)  &
                            - turb % v_mean(c) * turb % t_mean(c)
          wt_p(i) = wt_p(i) + turb % wt_res(c)  &
                            - turb % w_mean(c) * turb % t_mean(c)

          ut_mod(i) = ut_mod(i) + turb % ut_mean(c)
          vt_mod(i) = vt_mod(i) + turb % vt_mean(c)
          wt_mod(i) = wt_mod(i) + turb % wt_mean(c)

        end if
        n_count(i) = n_count(i) + 1
      end if
    end do
  end do


  ! Average over all processors
  do pl=1, n_prob-1
    call Comm_Mod_Global_Sum_Int(n_count(pl))

    call Comm_Mod_Global_Sum_Real(wall_p(pl))

    call Comm_Mod_Global_Sum_Real(tz_p(pl))
    call Comm_Mod_Global_Sum_Real(ti_p(pl))
    call Comm_Mod_Global_Sum_Real(w_p(pl))

    call Comm_Mod_Global_Sum_Real(uu_p(pl))
    call Comm_Mod_Global_Sum_Real(vv_p(pl))
    call Comm_Mod_Global_Sum_Real(ww_p(pl))
    call Comm_Mod_Global_Sum_Real(uw_p(pl))
    call Comm_Mod_Global_Sum_Real(kin_p(pl))
    call Comm_Mod_Global_Sum_Real(kin_mod_p(pl))

    count =  count + n_count(pl)

    if(Flow % heat_transfer) then
      call Comm_Mod_Global_Sum_Real(t_p(pl))
      call Comm_Mod_Global_Sum_Real(t2_p(pl))
      call Comm_Mod_Global_Sum_Real(t2_mod_p(pl))
      call Comm_Mod_Global_Sum_Real(ut_p(pl))
      call Comm_Mod_Global_Sum_Real(vt_p(pl))
      call Comm_Mod_Global_Sum_Real(wt_p(pl))
      call Comm_Mod_Global_Sum_Real(ut_mod(pl))
      call Comm_Mod_Global_Sum_Real(vt_mod(pl))
      call Comm_Mod_Global_Sum_Real(wt_mod(pl))
    end if
  end do

  call Comm_Mod_Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count(i)
      tz_p  (i) = tz_p (i)  / n_count(i)
      ti_p  (i) = ti_p (i)  / n_count(i)
      w_p   (i) = w_p  (i)  / n_count(i)

      uu_p(i) = uu_p(i)   / n_count(i)
      vv_p(i) = vv_p(i)   / n_count(i)
      ww_p(i) = ww_p(i)   / n_count(i)
      uw_p(i) = uw_p(i)   / n_count(i)
      kin_p(i) = kin_p(i) / n_count(i)
      kin_mod_p(i) = kin_mod_p(i) / n_count(i)

      if(Flow % heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        t2_p(i) = t2_p(i) / n_count(i)
        t2_mod_p(i) = t2_mod_p(i) / n_count(i)
        ut_p(i) = ut_p(i) / n_count(i)
        vt_p(i) = vt_p(i) / n_count(i)
        wt_p(i) = wt_p(i) / n_count(i)
        ut_mod(i) = ut_mod(i) / n_count(i)
        vt_mod(i) = vt_mod(i) / n_count(i)
        wt_mod(i) = wt_mod(i) / n_count(i)
      end if
    end if
  end do

  open(3, file = res_name)

  if(Flow % heat_transfer) then
    if(this_proc < 2) then
      write(*,'(a1,(a12, f12.6))')'#', ' Nu1 number = ',  &
               -tz_p(1) / (t_hot - t_cold)
      write(*,'(a1,(a12, f12.6))')'#', ' Nu2 number = ',  &
               (0.5/15.0)*(t_hot-ti_p(1)+ti_p(n_prob-1)-t_cold) / wall_p(1)
    end if
  end if

  write(i,'(a1,2x,a60)') '#', ' z,'                    //  &  !  1
                              ' u,'                    //  &  !  2
                              ' uu, vv, ww, uw'        //  &  !  3 -  6
                              ' kin'                   //  &  !  7
                              ' t, ut, vt, wt'                !  8 - 11

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(3,'(12es15.5e3)') wall_p(i),                     &  !  1
                              (ti_p(i) - t_cold)/t_diff,     &  !  2
                              (t_p(i)  - t_cold)/t_diff,     &  !  3
                              kin_p(i),                      &  !  4
                              kin_mod_p(i),                  &  !  5
                              (kin_p(i) + kin_mod_p(i)),     &  !  6
                              t2_p(i) / (t_diff*t_diff),     &  !  7
                              t2_mod_p(i) / (t_diff*t_diff), &  !  8 
                              (t2_p(i)+t2_mod_p(i))          &
                              / (t_diff*t_diff),             &  !  9
                              wt_p(i),                       &  ! 10
                              wt_mod(i),                     &  ! 11
                              wt_p(i) + wt_mod(i)               ! 12
    end if
  end do

  close(3)

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
  if(Flow % heat_transfer) then
    deallocate(t_p)
    deallocate(t2_p)
    deallocate(t2_mod_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  print *, '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
