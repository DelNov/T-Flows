!==============================================================================!
   subroutine User_Mod_Save_Swarm(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results for paerticles in homogeneous directions.             !
!                                                                              !
!   The results are then writen in files swarm_name_res.dat and                !
!   swarm_name_res_plus.dat                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: Flow
  type(Turb_Type),       target :: Turb
  type(Vof_Type),        target :: Vof
  type(Swarm_Type),      target :: Swarm
  integer, optional, intent(in) :: domain    ! current domain
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: n_prob, pl, c, i, n_points
  integer                  :: fu1, fu2, nb, nc
  character(SL)            :: coord_name
  character(SL)            :: swarm_result_name, swarm_result_name_plus
  real, allocatable        :: z_p(:), u_p(:), v_p(:), w_p(:),                &
                              uw_pp(:), uu_pp(:), vv_pp(:), ww_pp(:),        &
                              u_pp(:), v_pp(:), w_pp(:), vw_pp(:), uv_pp(:), &
                              ind(:), wall_p(:), y_plus_p(:), vis_t_p(:)
  integer, allocatable     :: n_p(:), n_count1(:), n_count2(:), n_states(:)
  real                     :: ubulk, re, cf_dean, cf, u_tau_p
  real                     :: density_const, visc_const
  logical                  :: there
!==============================================================================!

  ! Take aliases for the Flow 
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  ! constant fluid properties for single phase
  visc_const      = maxval(Flow % viscosity(:))
  density_const   = maxval(Flow % density(:))

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')

  if(.not. Swarm % statistics) return

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
  n_points = 0

  open(9, file=coord_name)

  ! Read the number of searching intervals
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(9,*) ind(pl), z_p(pl)
  end do
  close(9)

  ! Primary Flow arrays:
  allocate(n_p   (n_prob));  n_p    = 0
  allocate(wall_p(n_prob));  wall_p = 0.0
  allocate(u_p   (n_prob));  u_p    = 0.0
  allocate(v_p   (n_prob));  v_p    = 0.0
  allocate(w_p   (n_prob));  w_p    = 0.0

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Discrete phase arrays:
  allocate(u_pp (-nb:nc));  u_pp  = 0.0
  allocate(v_pp (-nb:nc));  v_pp  = 0.0
  allocate(w_pp (-nb:nc));  w_pp  = 0.0
  allocate(uu_pp(-nb:nc));  uu_pp = 0.0
  allocate(vv_pp(-nb:nc));  vv_pp = 0.0
  allocate(ww_pp(-nb:nc));  ww_pp = 0.0
  allocate(uw_pp(-nb:nc));  uw_pp = 0.0

  allocate(vis_t_p (-nb:nc));  vis_t_p  = 0.0
  allocate(y_plus_p(-nb:nc));  y_plus_p = 0.0

  allocate(n_states(n_prob));  n_states = 0
  allocate(n_count2(n_prob));  n_count2 = 0

  !----------------------!
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

        ! 2nd-moment-central-stationary statistics for swarm
        uu_pp(i) = uu_pp(i) + Swarm % uu(c)  &
                            - Swarm % u_mean(c) * Swarm % u_mean(c)
        vv_pp(i) = vv_pp(i) + Swarm % vv(c)  &
                            - Swarm % v_mean(c) * Swarm % v_mean(c)
        ww_pp(i) = ww_pp(i) + Swarm % ww(c)  &
                            - Swarm % w_mean(c) * Swarm % w_mean(c)
        uw_pp(i) = uw_pp(i) + Swarm % uw(c)  &
                            - Swarm % u_mean(c) * Swarm % w_mean(c)

        ! Averaging over the number of cells in this bin  
        n_states(i) = n_states(i) + 1
        !n_states(i) = n_states(i) + Swarm % n_states(c)
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

        vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)

        n_count2(i) = n_count2(i) + 1  ! counter 2 for the flowfield
      end if

    end do
  end do

  ! Average over all processors
  do pl = 1, n_prob-1

    ! Carrier Flow
    call Global % Sum_Int(n_count2(pl))
    call Global % Sum_Real(wall_p(pl))
    call Global % Sum_Real(u_p(pl))
    call Global % Sum_Real(v_p(pl))
    call Global % Sum_Real(w_p(pl))

    call Global % Sum_Real(vis_t_p (pl))
    call Global % Sum_Real(y_plus_p(pl))

    ! Swarm
    call Global % Sum_Int(n_states(pl))
    call Global % Sum_Real(u_pp(pl))
    call Global % Sum_Real(v_pp(pl))
    call Global % Sum_Real(w_pp(pl))
    call Global % Sum_Real(uu_pp(pl))
    call Global % Sum_Real(vv_pp(pl))
    call Global % Sum_Real(ww_pp(pl))
    call Global % Sum_Real(uw_pp(pl))
  end do

  call Global % Wait

  do i = 1, n_prob-1

    ! Background Flow
    if(n_count2(i) .ne. 0) then
      wall_p(i) = wall_p(i) / n_count2(i)
      u_p   (i) = u_p   (i) / n_count2(i)
      v_p   (i) = v_p   (i) / n_count2(i)
      w_p   (i) = w_p   (i) / n_count2(i)

      vis_t_p (i) = vis_t_p (i) / n_count2(i)
      y_plus_p(i) = y_plus_p(i) / n_count2(i)
    end if

    ! Swarm
    if(n_states(i) .ne. 0) then
      u_pp(i)  = u_pp(i) / (1.*n_states(i))
      v_pp(i)  = v_pp(i) / (1.*n_states(i))
      w_pp(i)  = w_pp(i) / (1.*n_states(i))

      uu_pp(i) = uu_pp(i) / (1.*n_states(i))
      vv_pp(i) = vv_pp(i) / (1.*n_states(i))
      ww_pp(i) = ww_pp(i) / (1.*n_states(i))
      uw_pp(i) = uw_pp(i) / (1.*n_states(i))
    end if
  end do

  ! Creating files for swarm statistics
  call File % Set_Name(swarm_result_name, appendix='-swarm-res',            &
                       time_step=Time % Curr_Dt(), extension='.dat')
  call File % Open_For_Writing_Ascii(swarm_result_name, fu1)
  call File % Set_Name(swarm_result_name_plus, appendix='-swarm-res-plus',  &
                       time_step=Time % Curr_Dt(), extension='.dat')
  call File % Open_For_Writing_Ascii(swarm_result_name_plus, fu2)

  ! Calculating friction velocity and friction temperature (for the flow!)
  u_tau_p = sqrt( (visc_const*sqrt(u_p(1)**2 +   &
                                   v_p(1)**2 +   &
                                   w_p(1)**2) / wall_p(1)))

  if(u_tau_p .eq. 0.0) then
    if(First_Proc()) then
      write(*,*) '# Friction velocity is zero in Save_Swarm.f90!'
    end if
    return
  end if

  open(fu1, file = swarm_result_name)
  open(fu2, file = swarm_result_name_plus)

  re      = density_const * ubulk * 2.0 / visc_const
  cf_dean = 0.073*(re)**(-0.25)
  cf      = u_tau_p**2/(0.5*ubulk**2)

  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk 
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', density_const * ubulk * 2.0/visc_const
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', density_const*u_tau_p/visc_const
  write(fu1,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
  write(fu1,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p 

  write(fu1,'(a1,2x,a50)') '#',' z,'                    //  &  !  1
                               ' u,'                    //  &  !  2
                               ' uu, vv, ww, uw'        //  &  !  3 - 6
                               ' kin'                          !  7

  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk 
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', density_const * ubulk * 2.0 / visc_const
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', density_const*u_tau_p / visc_const
  write(fu2,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0*(u_tau_p / ubulk)**2
  write(fu2,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p

  write(fu2,'(a1,2x,a50)') '#',' 1) z,'                             //  &
                               ' 2) u,'                             //  &
                               ' 3) uu, 4) vv, 5) ww, 6) uw'        //  &
                               ' 7) kin'

  do i = 1, n_prob
    if(n_count2(i) .ne. 0) then
      write(fu1,'(7es15.5e3)')  wall_p(i),                        & !  1
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

    u_pp   (i) = u_pp(i) / (u_tau_p + TINY)
    v_pp   (i) = v_pp(i) / (u_tau_p + TINY)
    w_pp   (i) = w_pp(i) / (u_tau_p + TINY)

    uu_pp (i) = uu_pp (i) / (u_tau_p**2 + TINY)
    vv_pp (i) = vv_pp (i) / (u_tau_p**2 + TINY)
    ww_pp (i) = ww_pp (i) / (u_tau_p**2 + TINY)
    uw_pp (i) = uw_pp (i) / (u_tau_p**2 + TINY)
  end do

  ! Writing normalized statistics 
  do i = 1, n_prob-1
    if(n_count2(i) .ne. 0) then
      write(fu2,'(7es15.5e3)')  wall_p(i),                        & !  1
                                u_pp(i),                          & !  2
                                uu_pp(i),                         & !  3
                                vv_pp(i),                         & !  4
                                ww_pp(i),                         & !  5
                                uw_pp(i),                         & !  6
                                0.5*(uu_pp(i)+vv_pp(i)+ww_pp(i))    !  7
    end if
  end do

  close(fu1)
  close(fu2)

  if(First_Proc())  print *, '# Finished with User_Mod_Save_Swarm.f90.'

  end subroutine
