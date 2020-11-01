!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm, n,    &
                                       n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function computes parasitic current intensities around a droplet      !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_v1  => r_cell_01,   &
                      dist_ck => r_Cell_02,   &
                      dist_cn => r_Cell_02,   &
                      dist_n  => i_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  integer                       :: n_stat_t, n_stat_p
  real                          :: time  ! physical time
!--------------------------------[Locals]--------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: s, c, c1, c2, last_cell, fu, n_tot_cells, c_dist
  real                     :: pos_mcl, h_drop !position mcl, droplet height
  real                     :: vol_wall_bot, vol_symm !volume of boundaries
  ! RMS u, area of droplet, area outside droplet
  ! pressure inside and outside, area of droplet
  real                     :: u_rms, a_in, a_out, p_in, p_out, a_vof
  real                     :: u_res, u_max, p_max, p_min
  real                     :: min_vfrac, max_vfrac
  real                     :: L_dom, H_dom ! width and jeight of domain
  real                     :: epsloc, norm_res, maxcurv, mincurv, ngrd
  real                     :: x0(3), x1(3), x2(3)
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  vof  => mult % vof

  epsloc = epsilon(epsloc)

  u_rms = 0.0
  a_in  = 0.0
  a_out = 0.0
  p_in  = 0.0
  p_out = 0.0
  a_vof = 0.0
  mincurv =  HUGE
  maxcurv = -HUGE

  ! Find max and min vfractions, to limit pressure calculation:

  min_vfrac = minval(vof % n(1:grid % n_cells - grid % comm % n_buff_cells))
  max_vfrac = maxval(vof % n(1:grid % n_cells - grid % comm % n_buff_cells))

  call Comm_Mod_Global_Min_Real(min_vfrac)
  call Comm_Mod_Global_Max_Real(max_vfrac)

  n_tot_cells = grid % n_cells - grid % comm % n_buff_cells

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    u_res = sqrt( flow % u % n(c) ** 2       &
                + flow % v % n(c) ** 2       &
                + flow % w % n(c) ** 2)
    sum_v1(c) = u_res
    u_rms = u_rms + u_res ** 2.0

    if (abs(max_vfrac - vof % n(c)) < epsloc) then
      a_in = a_in + grid % vol(c)
      p_in = p_in + flow % p % n(c) * grid % vol(c)
    end if

    if (abs(vof % n(c)) < min_vfrac + epsloc) then
      a_out = a_out + grid % vol(c)
      p_out = p_out + flow % p % n(c) * grid % vol(c)
    end if

    a_vof = a_vof + grid % vol(c) * vof % n(c)

    ngrd = norm2((/vof % x(c), vof % y(c), vof % z(c)/))
    if(ngrd > epsloc .and. mult % curv(c) > epsloc) then
      mincurv= min(mincurv, mult % curv(c))
      maxcurv= max(maxcurv, mult % curv(c))
    end if

  end do

  call Comm_Mod_Global_Sum_Real(a_vof)

  call Comm_Mod_Global_Sum_Int(n_tot_cells)
  call Comm_Mod_Global_Sum_Real(u_rms)
  u_rms = sqrt(1.0 / real(n_tot_cells) * u_rms)

  u_max = maxval(sum_v1(1:grid % n_cells - grid % comm % n_buff_cells))
  call Comm_Mod_Global_Max_Real(u_max)

  p_max = maxval(flow % p % n(1:grid % n_cells - grid % comm % n_buff_cells))
  p_min = minval(flow % p % n(1:grid % n_cells - grid % comm % n_buff_cells))
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)
  call Comm_Mod_Global_Max_Real(maxcurv)
  call Comm_Mod_Global_Min_Real(mincurv)

  call Comm_Mod_Global_Sum_Real(p_in)
  call Comm_Mod_Global_Sum_Real(p_out)
  call Comm_Mod_Global_Sum_Real(a_in)
  call Comm_Mod_Global_Sum_Real(a_out)

  p_in = p_in / a_in
  p_out = p_out / a_out

  ! Write Results
  if (this_proc < 2) then
    call File_Mod_Append_File_For_Writing('spurious.dat', fu)
    write(fu,'(8(2X,E16.10E2))')  time, a_vof, u_rms, u_max, &
                                  p_max-p_min, p_in-p_out, mincurv, maxcurv
    close(fu)
  end if

  end subroutine
