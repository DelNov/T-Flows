!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for dam break                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer                  :: n_stat_t, n_stat_p
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Var_Type),  pointer     :: fun
  integer                      :: c, last_cell, fu
  real                         :: b_volume, surface, rise_velocity,  &
                                  circularity, c_position
  integer                      :: s, p, c1, c2, closest_face, i_probe,  &
                                  code, i_s, th
  integer, parameter           :: N_PROBE = 4
  real                         :: vf_up, vf_down, z_up, z_down, height
  real                         :: h_probe(N_PROBE)
  real                         :: tmp_h(100)
  real                         :: dist, glo_dist, min_dist
  real,    contiguous, pointer :: p_node(:), vof_node(:)
  integer, contiguous, pointer :: face_dipped(:)
!==============================================================================!

  call Work % Connect_Real_Node(p_node, vof_node)
  call Work % Connect_Int_Face (face_dipped)

  ! Take aliases
  Grid  => Flow % pnt_grid
  fun   => Vof % fun

  ! Height probes:
  call Grid % Exchange_Cells_Real(Vof % fun % n)
  call Flow % Interpolate_Cells_To_Nodes(Vof % fun % n, vof_node)
  h_probe = 0.0

  ! Find fun at probe nodes
  do i_probe = 1, N_PROBE
    if (allocated(probes(i_probe) % s_probe)) then
      do i_s = 1, size(probes(i_probe) % s_probe)
        s = probes(i_probe) % s_probe(i_s)
        call Interpolate_From_Nodes(Grid, Vof, vof_node, s, i_probe, i_s)
      end do
    end if
  end do

  ! Find height by lineal interpolation on each probe
  do i_probe = 1, N_PROBE
    if (allocated(probes(i_probe) % s_probe)) then
      th = 0
      tmp_h = -HUGE
      do i_s = 2, size(probes(i_probe) % s_probe)
        if (      (probes(i_probe) % s_vof(i_s)     >= 0.5     &
            .and.  probes(i_probe) % s_vof(i_s - 1) < 0.5)     &
            .or.  (probes(i_probe) % s_vof(i_s)     <= 0.5     &
            .and.  probes(i_probe) % s_vof(i_s - 1) > 0.5) ) then
          th = th + 1
          tmp_h(th) =  probes(i_probe) % s_coor(i_s,3)        &
                    - ( probes(i_probe) % s_coor(i_s,3)       &
                      - probes(i_probe) % s_coor(i_s-1,3) )   &
                    * ( probes(i_probe) % s_vof(i_s) - 0.5 )  &
                    / ( probes(i_probe) % s_vof(i_s)          &
                      - probes(i_probe) % s_vof(i_s-1) )
        end if
      end do
      if (th > 0) then
        h_probe(i_probe) = maxval(tmp_h(1:th))
      end if
    end if
  end do

  do i_probe = 1, N_PROBE
    call Global % Max_Real(h_probe(i_probe))
  end do

  ! Pressure probes:
  call Grid % Exchange_Cells_Real(Flow % p % n)

  call Flow % Interpolate_Cells_To_Nodes(Flow % p % n, p_node)

  do i_probe = 1, size(nod_probe)
    if (nod_probe(i_probe) .ne. -1) then
      p_probe(i_probe) = p_node(nod_probe(i_probe))
    else
      p_probe(i_probe) = -HUGE
    end if
  end do

  do i_probe = 1, size(nod_probe)
    call Global % Max_Real(p_probe(i_probe))
  end do

  !---------------------!
  !   Volume of liquid  !
  !---------------------!
  b_volume = 0.0

  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * fun % n(c)
  end do

  call Global % Sum_Real(b_volume)

  ! Write to file
  if (First_Proc()) then
    call File % Append_For_Writing_Ascii('probe-data.dat', fu)
    write(fu,'((e16.10e2),9(2x,e16.10e2),4(2x,e16.10e2))')          &
          Time % Get_Time(), b_volume, p_probe(1:size(nod_probe)),  &
          h_probe(1:N_PROBE)
    close(fu)
  end if

  call Work % Disconnect_Real_Node(p_node, vof_node)
  call Work % Disconnect_Int_Face (face_dipped)

  end subroutine
