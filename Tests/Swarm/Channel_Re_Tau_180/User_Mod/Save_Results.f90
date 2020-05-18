!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, swarm, ts) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  use Const_Mod                      ! constants
  use Comm_Mod                       ! parallel stuff
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type
  use Bulk_Mod,  only: Bulk_Type
  use Var_Mod,   only: Var_Type
  use File_Mod,  only: problem_name
  use Turb_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Swarm_Type),      target :: swarm
  type(Multiphase_Type), target :: mult
  integer                       :: ts
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: n_prob, pl, c, i, count, s, c1, c2, n_points
  integer                  :: fu1, fu2
  character(len=80)        :: coord_name, result_name, result_name_plus
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:), t_p(:),  &
                              ind(:),  wall_p(:), kin_p(:), eps_p(:),  &
                              uw_p(:), uu_p(:), vv_p(:), ww_p(:),      &
                              t2_p(:), ut_p(:), vt_p(:), wt_p(:)
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: t_wall, t_tau, d_wall, nu_mean, t_inf
  real                     :: ubulk, error, re, cf_dean, cf, pr, u_tau_p
  real                     :: visc_const, dens_const, capa_const, cond_const
  logical                  :: there
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Field_Mod_Alias_Energy  (flow, t)
  visc_const = maxval(flow % viscosity(:))
  dens_const = maxval(flow % density(:))
  capa_const = maxval(flow % capacity(:))
  cond_const = maxval(flow % conductivity(:))

  ! Set the name for coordinate file
  call File_Mod_Set_Name(coord_name, extension='.1d')

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

  do c = 1, grid % n_cells
    ubulk    = bulk % flux_x / (dens_const*bulk % area_x)
    t_wall   = 0.0
    nu_mean  = 0.0
    n_points = 0
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

  allocate(n_p   (n_prob));  n_p      = 0
  allocate(wall_p(n_prob));  wall_p   = 0.0
  allocate(u_p   (n_prob));  u_p      = 0.0
  allocate(v_p   (n_prob));  v_p      = 0.0
  allocate(w_p   (n_prob));  w_p      = 0.0
  allocate(uu_p  (n_prob));  uu_p     = 0.0
  allocate(vv_p  (n_prob));  vv_p     = 0.0
  allocate(ww_p  (n_prob));  ww_p     = 0.0
  allocate(uw_p  (n_prob));  uw_p     = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0

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
    end if
  end do

  ! Calculating friction velocity and friction temperature
    do c = 1, grid % n_cells
      u_tau_p = sqrt( (visc_const*sqrt(u_p(1)**2 +   &
                                       v_p(1)**2 +   &
                                       w_p(1)**2)    &
                                       / wall_p(1))  &
                                       / dens_const)
    end do

  if(u_tau_p .eq. 0.0) then
    if(this_proc < 2) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if
    return
  end if

  call File_Mod_Set_Name(result_name, time_step = ts,              & 
       appendix='-res', extension='.dat')
  call File_Mod_Open_File_For_Writing(result_name, fu1)
  call File_Mod_Set_Name(result_name_plus, time_step = ts,         & 
       appendix='-res-plus', extension='.dat')
  call File_Mod_Open_File_For_Writing(result_name_plus, fu2)

  open(fu1,file=result_name)
  open(fu2,file=result_name_plus)


    pr = visc_const * capa_const / cond_const
    re = dens_const * ubulk * 2.0 / visc_const
    cf_dean = 0.073*(re)**(-0.25)
    cf      = u_tau_p**2/(0.5*ubulk**2)
    error   = abs(cf_dean - cf)/cf_dean * 100.0

    write(fu1,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(fu1,'(a1,(a12,e12.6))')  &
    '#', 're       = ', dens_const * ubulk * 2.0/visc_const
    write(fu1,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', dens_const*u_tau_p/visc_const
    write(fu1,'(a1,(a12,e12.6))')  &
    '#', 'cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(fu1,'(a1,(a12,f12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(fu1,'(a1,(a12,f12.6,a2,a22))') & 
    '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'


    write(fu2,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(fu2,'(a1,(a12,e12.6))')  &
    '#', 're       = ', dens_const * ubulk * 2.0/visc_const
    write(fu2,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', dens_const*u_tau_p/visc_const
    write(fu2,'(a1,(a12,e12.6))')  &
    '#', 'cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(fu2,'(a1,(a12,f12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(fu2,'(a1,(a12,f12.6,a2,a22))') & 
    '#', 'Cf_error = ', error, ' %', 'Dean formula is used.'

     
    write(fu1,'(a1,2x,a50)') '#',  ' z,'                    //  &  !  1
                                   ' u,'                    //  &  !  2
                                   ' uu, vv, ww, uw'        //  &  !  3 -  6
                                   ' kin'                          !  7

    write(fu2,'(a1,2x,a50)') '#',  ' z,'                    //  &  !  1
                                   ' u,'                    //  &  !  2
                                   ' uu, vv, ww, uw'        //  &  !  3 -  6
                                   ' kin'                          !  7
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu1,'(7es15.5e3)')  wall_p(i),                       & !  1
                                u_p(i),                          & !  2
                                uu_p(i),                         & !  3
                                vv_p(i),                         & !  4
                                ww_p(i),                         & !  5
                                uw_p(i),                         & !  6
                                0.5*(uu_p(i)+vv_p(i)+ww_p(i))      !  7
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
  end do

    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu2,'(7es15.5e3)')  wall_p(i),                     & !  1
                                u_p(i),                          & !  2
                                uu_p(i),                         & !  3
                                vv_p(i),                         & !  4
                                ww_p(i),                         & !  5
                                uw_p(i),                         & !  6
                                0.5*(uu_p(i)+vv_p(i)+ww_p(i))      !  7
      end if
    end do

  close(fu1)
  close(fu2)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(uu_p)
  deallocate(vv_p)
  deallocate(ww_p)
  deallocate(uw_p)

  if(this_proc < 2)  print *, '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
