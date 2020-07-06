!==============================================================================!
  subroutine User_Mod_Save_Swarm(flow, turb, swarm, save_name) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results for paerticles in homogeneous directions.             !
!                                                                              !
!   The results are then writen in files swarm_name_res.dat and                !
!   swarm_name_res_plus.dat                                                    !
!------------------------------------------------------------------------------!
  use Const_Mod                      ! constants
  use Comm_Mod                       ! parallel stuff
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type, heat_transfer, heat_flux, heat, &
                       capacity, conductivity, heated_area 
  use Bulk_Mod,  only: Bulk_Type
  use Var_Mod,   only: Var_Type
  use Turb_Mod 
  use Swarm_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  type(Swarm_Type), target :: swarm
  type(Particle_Type), pointer :: part
  character(len=*)         :: save_name
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: n_prob, pl, c, i, count, s, c1, c2, n_points, k
  integer                  :: index_p, j, ip, counter, blabla, ip0, ip01
  integer                  :: ip1, ip2, counter_k, l, n_ss, kk, counter_kk 
  integer                  :: label, counter1, counter2
  character(len=80)        :: coord_name, res_name, res_name_plus
  character(len=80)        :: swarm_res_name, swarm_res_name_plus
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:), t_p(:),      &
                              ind(:),  wall_p(:), kin_p(:), eps_p(:),      &
                              uw_pp(:), uu_pp(:), vv_pp(:), ww_pp(:),      &
                              u_pp(:), v_pp(:), w_pp(:), memo(:),          &
                              memo_z(:), memo_u(:), memo_v(:), memo_w(:),  &
                              memo_uu(:), memo_vv(:), memo_ww(:),          &
                              z_pp(:), memo_uw(:)
  integer, allocatable     :: n_p(:), n_count1(:), n_count2(:), store(:)
  integer, allocatable     :: n_states(:)
  real                     :: t_wall, t_tau, d_wall, nu_mean, t_inf
  real                     :: ubulk, re, cf_dean, cf, pr, u_tau_p, temp
  real                     :: temp_p1, temp_pp, temp_p2, slice, sslice_l 
  real                     :: sslice_u, sslice_diff, temp_kk
  real                     :: density_const, visc_const
  logical                  :: there, flag
!==============================================================================!


  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Field_Mod_Alias_Energy  (flow, t)

  ! constant fluid properties for single phase
  visc_const      = maxval(viscosity(:))
  density_const   = maxval(density(:))

  ! Set the name for coordinate file
  !call File_Mod_Set_Name(0, coord_name, ".1d")
  call File_Mod_Set_Name(coord_name, extension='.1d')

  call File_Mod_Set_Name(swarm_res_name, extension='-swarm-res.dat')
  call File_Mod_Set_Name(swarm_res_name_plus, extension='-swarm-res-plus.dat')

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

  ubulk    = bulk % flux_x / (density_const*bulk % area_x)
  t_wall   = 0.0
  nu_mean  = 0.0
  n_points = 0
  index_p  = 0
  temp     = 0.0
  temp_p1  = 0.0
  temp_p2  = 10.0**10
  temp_pp  = 10.0**10
  n_ss     = 4

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

  ! Primary flow arrays:
  allocate(n_p    (swarm % n_particles));   n_p      = 0
  allocate(wall_p (swarm % n_particles));   wall_p   = 0.0
  allocate(u_p    (swarm % n_particles));   u_p      = 0.0
  allocate(v_p    (swarm % n_particles));   v_p      = 0.0
  allocate(w_p    (swarm % n_particles));   w_p      = 0.0

  ! Discrete phase arrays:
  allocate(z_pp      (swarm % n_particles));   z_pp        = 0.0
  allocate(u_pp      (swarm % n_particles));   u_pp        = 0.0
  allocate(v_pp      (swarm % n_particles));   v_pp        = 0.0
  allocate(w_pp      (swarm % n_particles));   w_pp        = 0.0
  allocate(uu_pp     (swarm % n_particles));   uu_pp       = 0.0
  allocate(vv_pp     (swarm % n_particles));   vv_pp       = 0.0
  allocate(ww_pp     (swarm % n_particles));   ww_pp       = 0.0
  allocate(uw_pp     (swarm % n_particles));   uw_pp       = 0.0
  allocate(store     (swarm % n_particles));   store       = 0
  ! What if n_states>n_particles?!
  allocate(n_states  (swarm % n_particles));   n_states    = 0

  ! saving particle index history before arrangement
  allocate(memo_z    (swarm % n_particles));   memo_z      = 0.0
  allocate(memo_u    (swarm % n_particles));   memo_u      = 0.0
  allocate(memo_v    (swarm % n_particles));   memo_v      = 0.0
  allocate(memo_w    (swarm % n_particles));   memo_w      = 0.0
  allocate(memo_uu   (swarm % n_particles));   memo_uu     = 0.0
  allocate(memo_vv   (swarm % n_particles));   memo_vv     = 0.0
  allocate(memo_ww   (swarm % n_particles));   memo_ww     = 0.0
  allocate(memo_uw   (swarm % n_particles));   memo_uw     = 0.0

  allocate(n_count1(swarm % n_particles)); n_count1 = 0
  allocate(n_count2(swarm % n_particles)); n_count2 = 0
  count = 0

  !--------------- ------!
  !   Swarm statistics   !
  !----------------------!
  do i = 1, n_prob-1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells 
      if(grid % zc(c) > (z_p(i)) .and.  &
         grid % zc(c) < (z_p(i+1))) then

      ! Mean velocities
      u_pp(i) = u_pp(i) + swarm % u_mean(c)
      v_pp(i) = v_pp(i) + swarm % v_mean(c)
      w_pp(i) = w_pp(i) + swarm % w_mean(c)   

      ! 2nd-moment-central-stationary statistics for swarm
      uu_pp(i)    = uu_pp(i) + swarm % uu(c)  &
                  - swarm % u_mean(c) * swarm % u_mean(c)
      vv_pp(i)    = vv_pp(i) + swarm % vv(c)  &
                  - swarm % v_mean(c) * swarm % v_mean(c)
      ww_pp(i)    = ww_pp(i) + swarm % ww(c)  &
                  - swarm % w_mean(c) * swarm % w_mean(c)
      uw_pp(i)    = uw_pp(i) + swarm % uw(c)  &
                  - swarm % u_mean(c) * swarm % w_mean(c)

      ! Number of realizations
      n_states(i) = n_states(i) + swarm % n_states(c)   
      end if
    end do
  end do

  !--------------------------!
  !   Flowfield statistics   !
  !--------------------------!
  do i = 1, n_prob-1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells 
      if(grid % zc(c) > (z_p(i)) .and.  &
         grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + grid % wall_dist(c)

        ! Mean velocities
        u_p(i) = u_p(i) + turb % u_mean(c)
        v_p(i) = v_p(i) + turb % v_mean(c)
        w_p(i) = w_p(i) + turb % w_mean(c)

        n_count2(i) = n_count2(i) + 1  ! counter 2 for the flowfield 
      end if
    end do
  end do

  ! Average over all processors 
  do pl=1, n_prob-1

    call Comm_Mod_Global_Sum_Int(n_count2(pl))
    call Comm_Mod_Global_Sum_Real(wall_p(pl))
    call Comm_Mod_Global_Sum_Real(u_p(pl))
    call Comm_Mod_Global_Sum_Real(v_p(pl))
    call Comm_Mod_Global_Sum_Real(w_p(pl))

    counter2 =  counter2 + n_count2(pl)
  end do

  call Comm_Mod_Wait

  do i = 1, n_prob-1
    ! Background flow
    if(n_count2(i) .ne. 0) then
      wall_p(i)  = wall_p(i) / n_count2(i)
      u_p   (i)  = u_p   (i) / n_count2(i)
      v_p   (i)  = v_p   (i) / n_count2(i)
      w_p   (i)  = w_p   (i) / n_count2(i)
    end if 
   
    ! Particle-related 
    ! to be devided by n_states instead
    if(n_count1(i) .ne. 0) then

      u_pp(i) = u_pp(i) / n_states(i)
      v_pp(i) = v_pp(i) / n_states(i)
      w_pp(i) = w_pp(i) / n_states(i)

      uu_pp(i) = uu_pp(i) / n_states(i)
      vv_pp(i) = vv_pp(i) / n_states(i)
      ww_pp(i) = ww_pp(i) / n_states(i)
      uw_pp(i) = uw_pp(i) / n_states(i)

    end if
  end do



  ! Calculating friction velocity and friction temperature  (For the flow!)
    u_tau_p = sqrt( (visc_const*sqrt(u_p(1)**2 +   &
                                    v_p(1)**2 +   &
                                    w_p(1)**2)/ wall_p(1)))
  open(4, file = swarm_res_name_plus)

  do i = 3, 4
    pr = visc_const * capacity / conductivity
    re = density_const * ubulk * 2.0/visc_const
    cf_dean = 0.073*(re)**(-0.25)
    cf      = u_tau_p**2/(0.5*ubulk**2)
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'ubulk    = ', ubulk 
    write(i,'(a1,(a12,e12.6))')  &
    '#', 're       = ', density_const * ubulk * 2.0/visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'Re_tau   = ', density_const*u_tau_p/visc_const
    write(i,'(a1,(a12,e12.6))')  &
    '#', 'cf       = ', 2.0*(u_tau_p/ubulk)**2
    write(i,'(a1,(a12,f12.6))')  &
    '#', 'Utau     = ', u_tau_p 


      write(i,'(a1,2x,a50)') '#',  ' z,'                    //  &  !  1
                                   ' u,'                    //  &  !  2
                                   ' uu, vv, ww, uw'        //  &  !  3 -  6
                                   ' kin'                          !  7
  end do

    do i = 1, n_prob
      if(n_count2(i) .ne. 0) then   ! 16:45 22nd Nov
!      if(n_count1(i) .ne. 0) then   ! count index...
        write(3,'(7es15.5e3)')  wall_p(i),                        & !  1
!        write(3,'(7es15.5e3)')  z_pp(i),                        & !  1
                                u_pp(i),                          & !  2
                                uu_pp(i),                         & !  3
                                vv_pp(i),                         & !  4
                                ww_pp(i),                         & !  5
                                uw_pp(i),                         & !  6
                                0.5*(uu_pp(i)+vv_pp(i)+ww_pp(i))    !  7
      end if
    end do

  do i = 1, n_prob-1
    wall_p(i) = density_const * wall_p(i)*u_tau_p/visc_const

    z_pp   (i) = z_pp(i) / u_tau_p

    u_pp   (i) = u_pp(i) / u_tau_p
    v_pp   (i) = v_pp(i) / u_tau_p
    w_pp   (i) = w_pp(i) / u_tau_p

    uu_pp (i) = uu_pp (i) / (u_tau_p**2)
    vv_pp (i) = vv_pp (i) / (u_tau_p**2)
    ww_pp (i) = ww_pp (i) / (u_tau_p**2)
    uw_pp (i) = uw_pp (i) / (u_tau_p**2)
  end do

    !do i = 1, n_prob
    do i = 1, n_prob-1
      if(n_count2(i) .ne. 0) then
!      if(n_count1(i) .ne. 0) then
        write(4,'(7es15.5e3)')  wall_p(i),                          & !  1
!        write(4,'(7es15.5e3)')  z_pp(i),                          & !  1
                                u_pp(i),                            & !  2
                                uu_pp(i),                           & !  3
                                vv_pp(i),                           & !  4
                                ww_pp(i),                           & !  5
                                uw_pp(i),                           & !  6
                                0.5*(uu_pp(i)+vv_pp(i)+ww_pp(i))      !  7
      end if
    end do

  close(3)
  close(4)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(u_pp)
  deallocate(v_pp)
  deallocate(w_pp)
  deallocate(uu_pp)
  deallocate(vv_pp)
  deallocate(ww_pp)
  deallocate(uw_pp)

  deallocate(z_pp)
  deallocate(memo_z)
  deallocate(memo_u)
  deallocate(memo_v)
  deallocate(memo_w)
  deallocate(memo_uu)
  deallocate(memo_vv)
  deallocate(memo_ww)
  deallocate(memo_uw)
  deallocate(store)


  if(this_proc < 2)  print *, '# Finished with User_Mod_Save_Swarm.f90.'

  end subroutine
