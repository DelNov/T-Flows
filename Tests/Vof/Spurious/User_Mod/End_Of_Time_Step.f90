!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, turb, Vof, swarm, n,    &
                                       n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function computes position of contact line and high of droplet        !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_v1  => r_cell_01,   &
                      dist_ck => r_Cell_02,   &
                      dist_cn => r_Cell_02,   &
                      dist_n  => i_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n         ! time step
  integer, intent(in)      :: n_stat_t
  integer, intent(in)      :: n_stat_p
  real,    intent(in)      :: time      ! physical time
!--------------------------------[Locals]--------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: fun
  integer                  :: s, c, c1, c2, last_cell, fu, n_tot_cells, c_dist
  real                     :: pos_mcl, h_drop !position mcl, droplet height
  real                     :: vol_wall_bot, vol_symm !volume of boundaries
  ! RMS u, area of droplet, area outside droplet
  ! pressure inside and outside, area of droplet
  real                     :: u_rms, a_in, a_out, p_in, p_out, a_vof
  real                     :: u_res, u_max, p_max, p_min
  real                     :: min_vfrac, max_vfrac
  real                     :: L_dom, H_dom ! width and jeight of domain
  real                     :: epsloc, norm_res
  real                     :: x0(3), x1(3), x2(3)
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  fun  => Vof % fun

  epsloc = epsilon(epsloc)

  u_rms = 0.0
  a_in  = 0.0
  a_out = 0.0
  p_in  = 0.0
  p_out = 0.0
  a_vof = 0.0

  ! Find max and min vfractions, to limit pressure calculation:

  min_vfrac = minval(fun % n(1:grid % n_cells - grid % comm % n_buff_cells))
  max_vfrac = maxval(fun % n(1:grid % n_cells - grid % comm % n_buff_cells))

  call Comm_Mod_Global_Min_Real(min_vfrac)
  call Comm_Mod_Global_Max_Real(max_vfrac)

  n_tot_cells = grid % n_cells - grid % comm % n_buff_cells

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    u_res = sqrt( Flow % u % n(c) ** 2       &
                + Flow % v % n(c) ** 2       &
                + Flow % w % n(c) ** 2)
    sum_v1(c) = u_res
    u_rms = u_rms + u_res ** 2.0

    if (abs(max_vfrac - fun % n(c)) < epsloc) then
      a_in = a_in + grid % vol(c)
      p_in = p_in + Flow % p % n(c) * grid % vol(c)
    end if

    if (abs(fun % n(c)) < min_vfrac + epsloc) then
      a_out = a_out + grid % vol(c)
      p_out = p_out + Flow % p % n(c) * grid % vol(c)
    end if

    a_vof = a_vof + grid % vol(c) * fun % n(c)

  end do

  call Comm_Mod_Global_Sum_Real(a_vof)

  call Comm_Mod_Global_Sum_Int(n_tot_cells)
  call Comm_Mod_Global_Sum_Real(u_rms)
  u_rms = sqrt(1.0 / real(n_tot_cells) * u_rms)

  u_max = maxval(sum_v1(1:grid % n_cells - grid % comm % n_buff_cells))
  call Comm_Mod_Global_Max_Real(u_max)

  p_max = maxval(Flow % p % n(1:grid % n_cells - grid % comm % n_buff_cells))
  p_min = minval(Flow % p % n(1:grid % n_cells - grid % comm % n_buff_cells))
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  call Comm_Mod_Global_Sum_Real(p_in)
  call Comm_Mod_Global_Sum_Real(p_out)
  call Comm_Mod_Global_Sum_Real(a_in)
  call Comm_Mod_Global_Sum_Real(a_out)

  p_in = p_in / a_in
  p_out = p_out / a_out

  ! Write Results
  if (this_proc < 2) then
    if(n .eq. 1) then
      call File_Mod_Delete('spurious.dat')
    end if
    call File_Mod_Append_File_For_Writing('spurious.dat', fu)
    write(fu,'(6(2X,E16.10E2))')  time, a_vof, u_rms, u_max, &
                                  p_max-p_min, p_in-p_out
    close(fu)
  end if

  !! Calculate distance to center of cylinder
  !! For normals for curvature
  !if(this_proc <= 2) then
  !  dist_n = 0
  !  do c = 1, grid % n_cells
  !    norm_res = norm2((/fun % x(c), fun % y(c), fun % z(c)/))
  !    if(norm_res > epsloc) then
  !      x0 = (/1.0, grid % yc(c), 1.0/)
  !      x1 = (/grid % xc(c), grid % yc(c), grid % zc(c)/)
  !      x2 = x1 + 10.0 * (/Vof % fc_x(c), Vof % fc_y(c), Vof % fc_z(c)/) &
  !                       / norm2((/Vof % fc_x(c),   &
  !                                 Vof % fc_y(c),   &
  !                                 Vof % fc_z(c)/))
  !      dist_ck(c) = norm2( Math_Mod_Cross_Product(x0 - x1, x0 - x2)    &
  !                        / norm2(x2-x1) )

  !      dist_n(c) = 1
  !    end if
  !  end do

  !  call File_Mod_Open_File_For_Writing('dist_centk.dat', fu)

  !  do c = 1, grid % n_cells
  !    if(dist_n(c) == 1) then
  !      write(fu,'(4(2X,E16.10E2))') grid % xc(c), grid % yc(c), grid % zc(c), &
  !                                   dist_ck(c)
  !    end if
  !  end do

  !  close(fu)

  !! Calculate distance to center of cylinder
  !! For normals for normals
  !  dist_n = 0
  !  do c = 1, grid % n_cells
  !    norm_res = norm2((/fun % x(c), fun % y(c), fun % z(c)/))
  !    if(norm_res > epsloc) then
  !      x0 = (/1.0, grid % yc(c), 1.0/)
  !      x1 = (/grid % xc(c), grid % yc(c), grid % zc(c)/)

  !      x2 = x1 + 10.0 * (/Vof % fs_x(c), Vof % fs_y(c), Vof % fs_z(c)/) &
  !                       / norm2((/Vof % fs_x(c),   &
  !                                 Vof % fs_y(c),   &
  !                                 Vof % fs_z(c)/))

  !      dist_cn(c) = norm2( Math_Mod_Cross_Product(x0 - x1, x0 - x2)    &
  !                        / norm2(x2-x1) )
  !      dist_n(c) = 1
  !    end if
  !  end do

  !  call File_Mod_Open_File_For_Writing('dist_centn.dat', fu)

  !  do c = 1, grid % n_cells
  !    if(dist_n(c) == 1) then
  !      write(fu,'(4(2X,E16.10E2))') grid % xc(c), grid % yc(c), grid % zc(c), &
  !                                   dist_cn(c)
  !    end if
  !  end do

  !  close(fu)
  !end if

  end subroutine
