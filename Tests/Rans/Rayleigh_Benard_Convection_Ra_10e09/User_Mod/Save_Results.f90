!==============================================================================!
  subroutine User_Mod_Save_Results(flow, save_name) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  use Const_Mod                      ! constants
  use Comm_Mod                       ! parallel stuff
  use Grid_Mod,  only: Grid_Type
  use Grad_Mod
  use Field_Mod, only: Field_Type, heat_transfer, heat_flux,        &
                       density, viscosity, capacity, conductivity,  &
                       grav_x, grav_y, grav_z
  use Bulk_Mod,  only: Bulk_Type
  use Var_Mod,   only: Var_Type
  use Name_Mod,  only: problem_name
  use Rans_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  character(len=*)         :: save_name
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n_prob, pl, c, i, count, s, c1, c2, n_points
  character(len=80)   :: coord_name, res_name, res_name_plus
  character(len=80)   :: store_name
  real, allocatable   :: z_p(:), tz_p(:), ti_p(:), w_p(:), t_p(:),         &
                         y_plus_p(:),  ind(:),  wall_p(:), kin_p(:),       &
                         eps_p(:), kin_mod_p(:),                           &
                         uw_p(:), uu_p(:), vv_p(:), ww_p(:),               &
                         tt_p(:), ut_p(:), vt_p(:), wt_p(:), tt_mod_p(:),  &
                         ut_mod(:), vt_mod(:), wt_mod(:)
  integer, allocatable :: n_p(:), n_count(:)
  real                 :: t_wall, t_tau, d_wall, t_hot, t_cold, t_diff
  logical              :: there
!==============================================================================!

!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  ! Set some constants
  t_cold =  5.0
  t_hot  = 15.0
  t_diff =  t_hot - t_cold

  call Grad_Mod_Variable(t, .true.)
  
! Set the name for coordinate file
  call Name_File(0, coord_name, ".1d")

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(0, res_name,      "-res.dat")
  call Name_File(0, res_name_plus, "-res-plus.dat")

!  call Grad_Mod_For_Phi(grid, t % n, 3, phi_z, .true.)

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

    ! Restore the name and return
    problem_name = store_name
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
  if(heat_transfer) then
    allocate(t_p (n_prob));     t_p = 0.0
    allocate(tt_p(n_prob));     tt_p = 0.0
    allocate(tt_mod_p(n_prob)); tt_mod_p = 0.0
    allocate(ut_p(n_prob));     ut_p = 0.0
    allocate(vt_p(n_prob));     vt_p = 0.0
    allocate(wt_p(n_prob));     wt_p = 0.0
    allocate(ut_mod(n_prob));   ut_mod = 0.0
    allocate(vt_mod(n_prob));   vt_mod = 0.0
    allocate(wt_mod(n_prob));   wt_mod = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells 
      if(grid % zc(c) > (z_p(i)) .and.  &
         grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + grid % zc(c)
        tz_p  (i) = tz_p  (i) + t % z(c)
        ti_p  (i) = ti_p  (i) + t % n(c)
        w_p   (i) = w_p   (i) + w % mean(c)

        uu_p(i)   = uu_p(i) + uu % mean(c) - u % mean(c) * u % mean(c)
        vv_p(i)   = vv_p(i) + vv % mean(c) - v % mean(c) * v % mean(c)
        ww_p(i)   = ww_p(i) + ww % mean(c) - w % mean(c) * w % mean(c)
        uw_p(i)   = uw_p(i) + uw % mean(c) - u % mean(c) * w % mean(c)
        kin_p(i)  = kin_p(i) &
                    + 0.5*(uu % mean(c) - u % mean(c) * u % mean(c) &
                         + vv % mean(c) - v % mean(c) * v % mean(c) &
                         + ww % mean(c) - w % mean(c) * w % mean(c))
        kin_mod_p(i) = kin_mod_p(i) + kin % mean(c)  ! <--= we put kin % mean

        if(heat_transfer) then
          t_p(i)  = t_p(i)  + (t % mean(c) - t_cold)/t_diff
          tt_p(i) = tt_p(i) + tt % mean(c) - t % mean(c) * t % mean(c)
          tt_mod_p(i) = tt_mod_p(i) + t2 % mean(c)  ! <--= isn't it t2 % mean?
          ut_p(i) = ut_p(i) + ut % mean(c) - u % mean(c) * t % mean(c)
          vt_p(i) = vt_p(i) + vt % mean(c) - v % mean(c) * t % mean(c)
          wt_p(i) = wt_p(i) + wt % mean(c) - w % mean(c) * t % mean(c)
          ut_mod(i) = ut_mod(i) + ut % mean(c)  ! <--= isn't it ut % mean?
          vt_mod(i) = vt_mod(i) + vt % mean(c)  ! <--= isn't it vt % mean?
          wt_mod(i) = wt_mod(i) + wt % mean(c)  ! <--= isn't it wt % mean?
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

    if(heat_transfer) then
      call Comm_Mod_Global_Sum_Real(t_p(pl))
      call Comm_Mod_Global_Sum_Real(tt_p(pl))
      call Comm_Mod_Global_Sum_Real(tt_mod_p(pl))
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

      if(heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        tt_p(i) = tt_p(i) / n_count(i)
        tt_mod_p(i) = tt_mod_p(i) / n_count(i)
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

  if(heat_transfer) then
    if(this_proc < 2) then
      write(3,'(a1,(a12, f12.6))')'#', ' Nu number = ',  &
               tz_p(1) / (t_hot - t_cold)
      write(*,'(a1,(a12, f12.6))')'#', ' Nu number = ',  &
               0.05*(t_hot-ti_p(1)+ti_p(n_prob-1)-t_cold) / wall_p(1)
    end if
  end if

  write(i,'(a1,2x,a60)') '#', ' z,'                    //  &  !  1
                              ' u,'                    //  &  !  2
                              ' uu, vv, ww, uw'        //  &  !  3 -  6
                              ' kin'                   //  &  !  7
                              ' t, ut, vt, wt'                !  8 - 11

  do i = 1, n_prob-1
    t_p (i) = (t_p(i) - t_cold) / t_diff           ! t % n(c)
    tt_p(i) = tt_p(i) / (t_diff*t_diff)            ! ut % n(c)
    tt_mod_p(i) = tt_mod_p(i) / (t_diff*t_diff)    ! ut % n(c)
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(3,'(21es15.5e3)') wall_p(i),                  &  !  1
                              tz_p(i),                    &  !  2
                              (ti_p(i) - t_cold)/t_diff,  &  !  3
                              w_p(i),                     &  !  4
                              kin_p(i),                   &  !  5
                              kin_mod_p(i),               &  !  6
                              (kin_p(i) + kin_mod_p(i)),  &  !  7
                              uw_p(i),                    &  !  8
                              (t_p(i) - t_cold)/t_diff,   &  !  9
                              tt_p(i),                    &  ! 10
                              tt_mod_p(i),                &  ! 11
                              (tt_p(i)+tt_mod_p(i)),      &  ! 12
                              ut_p(i),                    &  ! 13
                              vt_p(i),                    &  ! 14
                              wt_p(i),                    &  ! 15
                              ut_mod(i),                  &  ! 16
                              vt_mod(i),                  &  ! 17
                              wt_mod(i),                  &  ! 18
                              ut_p(i) + ut_mod(i),        &  ! 19
                              vt_p(i) + vt_mod(i),        &  ! 20
                              wt_p(i) + wt_mod(i)            ! 21
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
  if(heat_transfer) then
    deallocate(t_p)
    deallocate(tt_p)
    deallocate(tt_mod_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  print *, '# Finished with User_Mod_Save_Results.f90.'

  ! Restore the name
  problem_name = store_name

  end subroutine
