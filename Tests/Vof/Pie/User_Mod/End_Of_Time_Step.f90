!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n_stat_t  ! 1st t.s. statistics turbulence
  integer, intent(in)      :: n_stat_p  ! 1st t.s. statistics particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell, fu
  real                     :: b_volume, surface, rise_vel_int
  real                     :: circularity, y_pos_cen
  real, save               :: rise_vel_cen, y_pos_cen_old
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  !----------------------------------------------!
  !   Bubble volume, surface and rise velocity   !
  !----------------------------------------------!
  b_volume = 0.0
  surface = 0.0
  y_pos_cen = 0.0
  rise_vel_int = 0.0
  call Flow % Grad_Variable(fun)

  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * fun % n(c)
    if (norm2((/fun % x(c),fun % y(c),fun % z(c)/)) > 1.0) then

      surface = surface + sqrt(fun % x(c) ** 2                    &
                             + fun % y(c) ** 2                    &
                             + fun % z(c) ** 2) * Grid % vol(c)
    end if
    y_pos_cen = y_pos_cen + Grid % yc(c) * fun % n(c) * Grid % vol(c)
    rise_vel_int = rise_vel_int + Flow % v % n(c) * fun % n(c) * Grid % vol(c)
  end do

  call Global % Sum_Real(b_volume)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(y_pos_cen)
  call Global % Sum_Real(rise_vel_int)
  y_pos_cen    = y_pos_cen    / b_volume
  rise_vel_int = rise_vel_int / b_volume
  rise_vel_cen = (y_pos_cen - y_pos_cen_old) / Flow % dt

  ! Just open the file benchmark.dat
  if(Time % Curr_Dt() .eq. 1) then
    call File % Open_For_Writing_Ascii('benchmark.dat', fu)
    close(fu)
  end if

  !-------------------!
  !   Write results   !
  !-------------------!
  if(Time % Curr_Dt() > 1) then

    if(First_Proc()) then
      print *, 'y_pos_cen        = ', y_pos_cen
      print *, 'y_pos_cen_old    = ', y_pos_cen_old
      print *, 'rise_vel_int (1) = ', rise_vel_int
      print *, 'rise_vel_cen (2) = ', rise_vel_cen

      ! Write to file
      call File % Append_For_Writing_Ascii('benchmark.dat', fu)

      ! With sphericity 3D:
      write(fu,'(6(2x,es16.10e2))') Time % Get_Time(), b_volume,             &
                                    PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)  &
                                    /surface,                                &
                                    y_pos_cen,                               &
                                    rise_vel_int,                            &
                                    rise_vel_cen
      close(fu)

    end if  ! First_Proc()

  end if  ! n > 1

  y_pos_cen_old = y_pos_cen

  end subroutine
