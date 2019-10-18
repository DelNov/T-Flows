!==============================================================================!
  subroutine Multiphase_Mod_Vof_Spurious_Post(mult, time_loc)
!------------------------------------------------------------------------------!
!   Computes quantities of interest                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_v1 => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: time_loc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: s, c, cellt
  real                      :: u_max, p_min, p_max, u_rms
  real                      :: a_in, a_out, a_tot, u_res, p_in, p_out, a_vof
  real                      :: min_x, max_z
!==============================================================================!

  ! Take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid

  min_x =  HUGE
  max_z = -HUGE

  u_rms = 0.0
  a_tot = 0.0
  a_in  = 0.0
  a_out = 0.0
  p_in  = 0.0
  p_out = 0.0
  a_vof = 0.0

  do c = 1, grid % n_cells
    u_res = sqrt( flow % u % n(c) ** 2.0       &
                + flow % v % n(c) ** 2.0       &
                + flow % w % n(c) ** 2.0)
    sum_v1(c) = u_res 
    u_rms = u_rms + u_res ** 2.0
    a_tot = a_tot + grid % vol(c)
    if (abs(1.0 - mult % vof % n(c)) < TINY) then
      a_in = a_in + grid % vol(c)
      p_in = p_in + flow % p % n(c) * grid % vol(c) 
    end if

    if (abs(mult % vof % n(c)) < TINY) then
      a_out = a_out + grid % vol(c)
      p_out = p_out + flow % p % n(c) * grid % vol(c) 
    end if

    a_tot = a_tot + grid % vol(c)
    a_vof = a_vof + grid % vol(c) * mult % vof % n(c)

    if ((grid % xc(c) <= min_x)) then
      if ((grid % zc(c) >= max_z)) then
       cellt = c
       min_x = grid % xc(c)
       max_z = grid % zc(c)

     end if
    end if
  end do

  u_rms = sqrt(1.0 / real(grid % n_cells) * u_rms)
  u_max = maxval(sum_v1(:))
  p_max = maxval(flow % p % n(:))
  p_min = minval(flow % p % n(:))
  p_in = p_in / a_in
  p_out = p_out / a_out

!  write(*,*) cellt,grid % xc(cellt),grid % yc(cellt),grid % zc(cellt)
!  write(*,*) cellt,flow % p % n(cellt),flow % u %n(cellt),flow % v %n(cellt),flow % w %n(cellt)
!  cellt =maxloc(sum_v1(:),dim = 1)
!  write(*,*) cellt,sum_v1(cellt),grid % xc(cellt),grid % yc(cellt),grid % zc(cellt)
  open(9, file = 'Spurious-post.dat',position='append')

    write(9,*) time_loc, a_vof, u_rms, u_max, p_max-p_min, p_in-p_out

  close(9)

  end subroutine
