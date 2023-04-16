!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!                                                                              !
!   In its current form, it only extracts the results in specified probes      !
!   prints values of mean velocity and mean Reynolds stresses and exits.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  integer                      :: n, n_probes, backup_dt
  real                         :: x_p, y_p     ! coordinates of the probe
  real                         :: backup_time
  character(80)                :: arg
  real,    contiguous, pointer :: d_probe(:), z_probe(:)
  real,    contiguous, pointer :: u_mean_n(:), uu_res_n(:), vv_res_n(:)
  integer, contiguous, pointer :: node_ind(:)
!==============================================================================!

  call Work % Connect_Real_Node(d_probe, z_probe, u_mean_n, uu_res_n, vv_res_n)
  call Work % Connect_Int_Node (node_ind)

  Grid => Flow % pnt_grid

  !--------------------------!
  !   Check the invocation   !
  !--------------------------!
  if(command_argument_count() .eq. 2) then

    ! Get x_p and y_p
    call get_command_argument(1, arg);  read(arg, *) x_p
    call get_command_argument(2, arg);  read(arg, *) y_p
  else
    print *, '# You failed to invoke the program properly.'
    print *, '# Correct invocation:'
    print *, './Process  x_probe  y_probe'
    stop
  end if

  !-----------------------------------------------!
  !   Find nodes closest to the specified probe   !
  !-----------------------------------------------!

  ! Store node indices, distances to the probe and nodes' z coordinates
  d_probe(:) = HUGE
  do n = 1, Grid % n_nodes
    node_ind(n) = n
    d_probe(n) = min(d_probe(n),  &
                     sqrt((Grid % xn(n)-x_p)**2 + (Grid % yn(n)-y_p)**2) )
    z_probe(n) = Grid % zn(n)
  end do

  ! Sort nodes by their distance from probe and z coordinate
  call Sort % Two_Real_Carry_Int(d_probe, z_probe, node_ind)

  ! Find number of probes (use the z coordinate for that)
  do n = 1, Grid % n_nodes-1
    if(z_probe(n+1) < z_probe(n)) then
      n_probes = n
      goto 1
    end if
  end do
1 continue

  !--------------------------------------------------!
  !   Interpolate cell-based to node based results   !
  !--------------------------------------------------!
  call Flow % Interpolate_Cells_To_Nodes(Turb % u_mean, u_mean_n)
  call Flow % Interpolate_Cells_To_Nodes(Turb % uu_res, uu_res_n)
  call Flow % Interpolate_Cells_To_Nodes(Turb % vv_res, vv_res_n)

  !---------------------------!
  !   Print the profile out   !
  !---------------------------!
  print *, '# Profile at (x,y) = ', x_p, y_p
  do n = 1, n_probes
    print *, z_probe(n), u_mean_n(node_ind(n)),                             &
                         uu_res_n(node_ind(n)) - u_mean_n(node_ind(n))**2,  &
                         vv_res_n(node_ind(n))
  end do

  call Work % Disconnect_Real_Node(d_probe,z_probe,u_mean_n,uu_res_n,vv_res_n)
  call Work % Disconnect_Int_Node (node_ind)

  call Global % End_Parallel
  stop

  end subroutine
