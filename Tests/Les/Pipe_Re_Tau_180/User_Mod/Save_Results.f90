!==============================================================================!
  subroutine User_Mod_Save_Results(grid, save_name) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  use Grid_Mod
  use Flow_Mod
  use Rans_Mod
  use Comm_Mod                       ! parallel stuff
  use Name_Mod,  only: problem_name
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: save_name
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n_prob, pl, c, i, count, s, c1, c2, n_points
  character(len=80)   :: coord_name, res_name, res_name_plus
  character(len=80)   :: store_name
  real,allocatable    :: z_p(:), u_p(:), v_p(:), w_p(:), t_p(:), y_plus_p(:),  &
                         ind(:),  wall_p(:), kin_p(:), eps_p(:),               &
                         uw_p(:), uu_p(:), vv_p(:), ww_p(:),                   &
                         tt_p(:), ut_p(:), vt_p(:), wt_p(:),                   &
  integer,allocatable :: n_p(:), n_count(:)
  real                :: t_wall, t_tau, d_wall, nu_max 
  real                :: ubulk, error, re, cf_dean, cf, pr, u_tau_p
  logical             :: there
!==============================================================================!

  ! Set the name for coordinate file
  call Name_File(0, coord_name, ".1d")

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(0, res_name,      "-res.dat")
  call Name_File(0, res_name_plus, "-res-plus.dat")

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file=coord_name, exist=there)
  if(.not. there) then
    if(this_proc < 2) then
      print *, '==============================================================='
      print *, 'In order to extract profiles and write them in ascii files'
      print *, 'the code has to read cell-faces coordinates '
      print *, 'in wall-normal direction in the ascii file ''case_name.1d.'''
      print *, 'The file format should be as follows:'
      print *, '10  ! number of cells + 1'
      print *, '1 0.0'
      print *, '2 0.1'
      print *, '3 0.2'
      print *, '... '
      print *, '==============================================================='
    end if

    ! Restore the name and return
    problem_name = store_name
    return
  end if

  ubulk = bulk % flux_x / (density*bulk % area_x)
  t_wall = 0.0
  nu_max = 0.0
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

  allocate(n_p   (n_prob));  n_p      = 0
  allocate(wall_p(n_prob));  wall_p   = 0.0
  allocate(u_p   (n_prob));  u_p      = 0.0
  allocate(v_p   (n_prob));  v_p      = 0.0
  allocate(w_p   (n_prob));  w_p      = 0.0
  allocate(uu_p  (n_prob));  uu_p     = 0.0
  allocate(vv_p  (n_prob));  vv_p     = 0.0
  allocate(ww_p  (n_prob));  ww_p     = 0.0
  allocate(uw_p  (n_prob));  uw_p     = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(heat_transfer) then
    allocate(t_p (n_prob));  t_p = 0.0
    allocate(tt_p(n_prob));  tt_p = 0.0
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
        u_p   (i) = u_p   (i) + u % mean(c)
        v_p   (i) = v_p   (i) + v % mean(c)
        w_p   (i) = w_p   (i) + w % mean(c)

        uu_p(i) = uu_p(i) + uu % mean(c) - u % mean(c) * u % mean(c)
        vv_p(i) = vv_p(i) + vv % mean(c) - v % mean(c) * v % mean(c)
        ww_p(i) = ww_p(i) + ww % mean(c) - w % mean(c) * w % mean(c)
        uw_p(i) = uw_p(i) + uw % mean(c) - u % mean(c) * w % mean(c)

        if(heat_transfer) then
          t_p(i)  = t_p(i)  + t % mean(c)
          tt_p(i) = tt_p(i) + tt % mean(c) - t % mean(c) * t % mean(c)
          ut_p(i) = ut_p(i) + ut % mean(c) - u % mean(c) * t % mean(c)
          vt_p(i) = vt_p(i) + vt % mean(c) - v % mean(c) * t % mean(c)
          wt_p(i) = wt_p(i) + wt % mean(c) - w % mean(c) * t % mean(c)
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

    call Comm_Mod_Global_Sum_Real(uu_p(pl))
    call Comm_Mod_Global_Sum_Real(vv_p(pl))
    call Comm_Mod_Global_Sum_Real(ww_p(pl))
    call Comm_Mod_Global_Sum_Real(uw_p(pl))

    count =  count + n_count(pl)

    if(heat_transfer) then
      call Comm_Mod_Global_Sum_Real(t_p(pl))
      call Comm_Mod_Global_Sum_Real(tt_p(pl))
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

      uu_p(i) = uu_p(i) / n_count(i)
      vv_p(i) = vv_p(i) / n_count(i)
      ww_p(i) = ww_p(i) / n_count(i)
      uw_p(i) = uw_p(i) / n_count(i)

      if(heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        tt_p(i) = tt_p(i) / n_count(i)
        ut_p(i) = ut_p(i) / n_count(i)
        vt_p(i) = vt_p(i) / n_count(i)
        wt_p(i) = wt_p(i) / n_count(i)
      end if
    end if
  end do

  ! Calculating friction velocity and friction temperature
    u_tau_p = sqrt( (viscosity*sqrt(u_p(1)**2 +   &
                                    v_p(1)**2 +   &
                                    w_p(1)**2)    &
                                    / wall_p(1))  &
                                    / density)
  if(u_tau_p .eq. 0.0) then
    if(this_proc < 2) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if

    ! Restore the name and return
    problem_name = store_name
    return
  end if

  if(heat_transfer) then 
    d_wall = 0.0 
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) > d_wall) then
        d_wall = grid % wall_dist(c)
        t_inf  = t % mean(c)
      end if
    end do

    call Comm_Mod_Wait

    if(heat_flux> 0.0) then
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

          t_wall = t_wall + t % mean(c2)
          nu_max = nu_max + t % q(c2)/(conductivity*(t % mean(c2) - t_inf))
          n_points = n_points + 1
        end if
      end if
    end do

    call Comm_Mod_Global_Sum_Real(t_wall)
    call Comm_Mod_Global_Sum_Real(nu_max)
    call Comm_Mod_Global_Sum_Int(n_points)

    call Comm_Mod_Wait

    t_wall = t_wall / n_points
    nu_max = nu_max / n_points
    t_tau  = heat_flux / (density * capacity * u_tau_p)
  end if

  open(3, file = res_name)
  open(4, file = res_name_plus)
  
  do i = 3, 4
    pr = viscosity * capacity / conductivity
    re = density * ubulk * 2.0/viscosity
    cf_dean = 0.073*(re)**(-0.25)
    cf      = u_tau_p**2/(0.5*ubulk**2)
    error   = abs(cf_dean - cf)/cf_dean * 100.0
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(i,'(a1,(a12,e12.6))')  &
    '#', 're       = ', density * ubulk * 2.0/viscosity
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', density*u_tau_p/viscosity
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(i,'(a1,(a12,f12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(i,'(a1,(a12,f12.6,a2,a22))') & 
    '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'
    if(heat_transfer) then
      write(i,'(a1,(a12, f12.6))')'#', 'Nu number =', nu_max 
      write(i,'(a1,(a12, f12.6,a2,a39))')'#', 'Nu_error  =',  &
            abs(0.023*0.5*re**0.8*pr**0.4 - nu_max)           &
            / (0.023*0.5*re**0.8*pr**0.4) * 100.0, ' %',      &
            'correlation of Dittus-Boelter is used.' 
    end if

    if(heat_transfer) then
      write(i,'(a1,2x,a60)') '#',  ' z,'                    //  &
                                   ' u,'                    //  &
                                   ' uu, vv, ww, uw'        //  &
                                   ' kin'                   //  &
                                   ' t, ut, vt, wt,'   
    else
      write(i,'(a1,2x,a50)') '#',  ' z,'                    //  &
                                   ' u,'                    //  &
                                   ' uu, vv, ww, uw'        //  &
                                   ' kin'  
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(12e15.7)') wall_p(i),                       &
                             u_p(i),                          &
                             uu_p(i),                         &
                             vv_p(i),                         &
                             ww_p(i),                         &
                             uw_p(i),                         &
                             0.5*(uu_p(i)+vv_p(i)+ww_p(i)),   &
                             t_p(i),                          &
                             tt_p(i),                         &
                             ut_p(i),                         &
                             vt_p(i),                         &
                             wt_p(i)
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(7e15.7)')  wall_p(i),                       &
                             u_p(i),                          &
                             uu_p(i),                         &
                             vv_p(i),                         &
                             ww_p(i),                         &
                             uw_p(i),                         &
                             0.5*(uu_p(i)+vv_p(i)+ww_p(i))
      end if
    end do
  end if

  do i = 1, n_prob-1
    wall_p(i) = density * wall_p(i)*u_tau_p/viscosity
    u_p   (i) = u_p(i) / u_tau_p
    v_p   (i) = v_p(i) / u_tau_p
    w_p   (i) = w_p(i) / u_tau_p

    kin_p(i) = kin_p(i) / u_tau_p**2                      ! kin%n(c)
    eps_p(i) = eps_p(i)*viscosity / (u_tau_p**4*density)  ! eps%n(c)
    uu_p (i) = uu_p (i) / (u_tau_p**2)
    vv_p (i) = vv_p (i) / (u_tau_p**2)
    ww_p (i) = ww_p (i) / (u_tau_p**2)
    uw_p (i) = uw_p (i) / (u_tau_p**2)

    if(heat_transfer) then
      t_p (i) = (t_wall - t_p(i)) / t_tau  ! t % n(c)
      tt_p(i) = tt_p(i) / (t_tau*t_tau)    ! ut % n(c)
      ut_p(i) = ut_p(i) / (u_tau_p*t_tau)  ! ut % n(c)
      vt_p(i) = vt_p(i) / (u_tau_p*t_tau)  ! vt % n(c)
      wt_p(i) = wt_p(i) / (u_tau_p*t_tau)  ! wt % n(c)
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(12e15.7)') wall_p(i),                       &
                             u_p(i),                          &
                             uu_p(i),                         &
                             vv_p(i),                         &
                             ww_p(i),                         &
                             uw_p(i),                         &
                             0.5*(uu_p(i)+vv_p(i)+ww_p(i)),   &
                             t_p(i),                          &
                             tt_p(i),                         &
                             ut_p(i),                         &
                             vt_p(i),                         &
                             wt_p(i)
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(7e15.7)')  wall_p(i),                       &
                             u_p(i),                          &
                             uu_p(i),                         &
                             vv_p(i),                         &
                             ww_p(i),                         &
                             uw_p(i),                         &
                             0.5*(uu_p(i)+vv_p(i)+ww_p(i))
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
  if(heat_transfer) then
    deallocate(t_p)
    deallocate(tt_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  ! Restore the name
  problem_name = store_name

  end subroutine
