!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results for paerticles in homogeneous directions.             !
!                                                                              !
!   The results are then writen in files swarm_name_res.dat and                !
!   swarm_name_res_plus.dat                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target  :: Flow
  type(Turb_Type),     target  :: Turb
  type(Vof_Type),      target  :: Vof
  type(Swarm_Type),    target  :: Swarm
  integer, optional            :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: n_prob, pl, c, i, count, s, c1, c2, n_points, k
  integer                  :: index_p, j, counter
  integer                  :: counter_k, l, n_ss, kk, counter_kk
  integer                  :: label, counter1, counter2
  character(SL)            :: coord_name, res_name, res_name_plus
  character(SL)            :: swarm_res_name, swarm_res_name_plus
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
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  ! Constant fluid properties for single phase
  visc_const    = maxval(Flow % viscosity(:))
  density_const = maxval(Flow % density(:))

  ! Set the name for coordinate file
  !call File % Set_Name(0, coord_name, ".1d")
  call File % Set_Name(coord_name, extension='.1d')

  call File % Set_Name(swarm_res_name, extension='-Swarm-res.dat')
  call File % Set_Name(swarm_res_name_plus, extension='-Swarm-res-plus.dat')

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

  ! Primary Flow arrays:
  allocate(n_p    (Swarm % n_particles));   n_p      = 0
  allocate(wall_p (Swarm % n_particles));   wall_p   = 0.0
  allocate(u_p    (Swarm % n_particles));   u_p      = 0.0
  allocate(v_p    (Swarm % n_particles));   v_p      = 0.0
  allocate(w_p    (Swarm % n_particles));   w_p      = 0.0

  ! Discrete phase arrays:
  allocate(z_pp      (Swarm % n_particles));   z_pp        = 0.0
  allocate(u_pp      (Swarm % n_particles));   u_pp        = 0.0
  allocate(v_pp      (Swarm % n_particles));   v_pp        = 0.0
  allocate(w_pp      (Swarm % n_particles));   w_pp        = 0.0
  allocate(uu_pp     (Swarm % n_particles));   uu_pp       = 0.0
  allocate(vv_pp     (Swarm % n_particles));   vv_pp       = 0.0
  allocate(ww_pp     (Swarm % n_particles));   ww_pp       = 0.0
  allocate(uw_pp     (Swarm % n_particles));   uw_pp       = 0.0
  allocate(store     (Swarm % n_particles));   store       = 0
  ! What if n_states>n_particles?!
  allocate(n_states  (Swarm % n_particles));   n_states    = 0

  ! saving particle index history before arrangement
  allocate(memo_z    (Swarm % n_particles));   memo_z      = 0.0
  allocate(memo_u    (Swarm % n_particles));   memo_u      = 0.0
  allocate(memo_v    (Swarm % n_particles));   memo_v      = 0.0
  allocate(memo_w    (Swarm % n_particles));   memo_w      = 0.0
  allocate(memo_uu   (Swarm % n_particles));   memo_uu     = 0.0
  allocate(memo_vv   (Swarm % n_particles));   memo_vv     = 0.0
  allocate(memo_ww   (Swarm % n_particles));   memo_ww     = 0.0
  allocate(memo_uw   (Swarm % n_particles));   memo_uw     = 0.0

  allocate(n_count1(Swarm % n_particles)); n_count1 = 0
  allocate(n_count2(Swarm % n_particles)); n_count2 = 0
  count = 0

  !--------------- ------!
  !   Swarm statistics   !
  !----------------------!
  do i = 1, n_prob-1
    do c = Cells_In_Domain()
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

      ! Mean velocities
      u_pp(i) = u_pp(i) + Swarm % u_mean(c)
      v_pp(i) = v_pp(i) + Swarm % v_mean(c)
      w_pp(i) = w_pp(i) + Swarm % w_mean(c)

      ! 2nd-moment-central-stationary statistics for Swarm
      uu_pp(i)    = uu_pp(i) + Swarm % uu(c)  &
                  - Swarm % u_mean(c) * Swarm % u_mean(c)
      vv_pp(i)    = vv_pp(i) + Swarm % vv(c)  &
                  - Swarm % v_mean(c) * Swarm % v_mean(c)
      ww_pp(i)    = ww_pp(i) + Swarm % ww(c)  &
                  - Swarm % w_mean(c) * Swarm % w_mean(c)
      uw_pp(i)    = uw_pp(i) + Swarm % uw(c)  &
                  - Swarm % u_mean(c) * Swarm % w_mean(c)

      ! Number of realizations
      n_states(i) = n_states(i) + Swarm % n_states(c)
      end if
    end do
  end do

  !--------------------------!
  !   Flowfield statistics   !
  !--------------------------!
  do i = 1, n_prob-1
    do c = Cells_In_Domain()
      if(Grid % zc(c) > (z_p(i)) .and.  &
         Grid % zc(c) < (z_p(i+1))) then

        wall_p(i) = wall_p(i) + Grid % wall_dist(c)

        ! Mean velocities
        u_p(i) = u_p(i) + Turb % u_mean(c)
        v_p(i) = v_p(i) + Turb % v_mean(c)
        w_p(i) = w_p(i) + Turb % w_mean(c)

        n_count2(i) = n_count2(i) + 1  ! counter 2 for the flowfield
      end if
    end do
  end do

  ! Average over all processors
  do pl=1, n_prob-1

    call Global % Sum_Int(n_count2(pl))
    call Global % Sum_Real(wall_p(pl))
    call Global % Sum_Real(u_p(pl))
    call Global % Sum_Real(v_p(pl))
    call Global % Sum_Real(w_p(pl))

    counter2 =  counter2 + n_count2(pl)
  end do

  call Global % Wait

  do i = 1, n_prob-1
    ! Background Flow
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
    pr = visc_const * maxval(Flow % capacity) / maxval(Flow % conductivity)
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
!   if(n_count1(i) .ne. 0) then   ! count index...
      write(3,'(7es15.5e3)')  wall_p(i),                        & !  1
!     write(3,'(7es15.5e3)')  z_pp(i),                          & !  1
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
!   if(n_count1(i) .ne. 0) then
      write(4,'(7es15.5e3)')  wall_p(i),                          & !  1
!     write(4,'(7es15.5e3)')  z_pp(i),                            & !  1
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


  if(First_Proc())  print *, '# Finished with User_Mod_Save_Swarm.f90.'

  end subroutine
