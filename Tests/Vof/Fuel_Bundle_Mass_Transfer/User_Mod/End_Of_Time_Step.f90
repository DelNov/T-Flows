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
  integer                  :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                  :: n_stat_p  ! 1st t.s. statistics particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell, fu
  real,parameter           :: xy_cent = 6.4e-3
  real,parameter           :: half_range = 1.5e-3
  real                     :: range_min, range_max
  real                     :: vol_l,vol_v, vol_vb, sum_vap, sum_cond
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  !-------------------------!
  !   volume of each phase  !
  !-------------------------!
  vol_l = 0.0   ! volume of liquid
  vol_v = 0.0   ! volume of vapor

  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    vol_l = vol_l + Grid % vol(c) * fun % n(c)
    vol_v = vol_v + Grid % vol(c) * (1.0-fun % n(c))
  end do

  call Global % Sum_Real(vol_l)
  call Global % Sum_Real(vol_v)

  !---------------------------!
  !   volume of vapor in box  !
  !---------------------------!
  vol_vb = 0.0  ! volume of vapor in box
  range_min = xy_cent - half_range
  range_max = xy_cent + half_range

  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    if (range_min < Grid % xc(c) .AND. Grid % xc(c) < range_max .AND.  &
        range_min < Grid % yc(c) .AND. Grid % yc(c) < range_max ) then
    vol_vb = vol_vb + Grid % vol(c) * (1.0-fun % n(c))
    endif
  end do

  call Global % Sum_Real(vol_vb)


  ! Write to file
  if (First_Proc()) then
    call File % Append_For_Writing_Ascii('volume.out', fu)

    ! With sphericity 3D
    write(fu,'(4(2X,E18.10E2))') Time % Get_Time(), vol_l, vol_v, vol_vb
    close(fu)
  end if

  !----------------------------------------!
  !  sum of vaporization and condensation  !
  !----------------------------------------!
  sum_vap=0.0    ! sum of vaporization [kg/s]
  sum_cond=0.0   ! sum of condensation [kg/s]
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    if (Vof % m_dot(c) > 0.0) then
      sum_vap = sum_vap + Vof % m_dot(c)
    else
      sum_cond = sum_cond + Vof % m_dot(c)
    endif
  end do

  call Global % Sum_Real(sum_vap)
  call Global % Sum_Real(sum_cond)

  ! Write to file
  if (First_Proc()) then
    call File % Append_For_Writing_Ascii('massTransfer.out', fu)

    ! With sphericity 3D
    write(fu,'(3E18.10E2)') Time % Get_Time(),sum_vap, sum_cond 
    close(fu)
  end if


  end subroutine
