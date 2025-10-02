!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: s, c1, c2
  real                     :: vol_in, area_in, vel_in
  real                     :: vel_max, vel_min
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  !-------------------------------------!
  !   Compute average inflow velocity   !
  !-------------------------------------!
  vol_in  = 0.0
  area_in = 0.0

  vel_max = -HUGE
  vel_min = +HUGE
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if( Grid % Bnd_Cond_Type(c2) .eq. PRESSURE ) then
        area_in = area_in + Grid % s(s)
        vol_in  = vol_in  + Grid % s(s) * Flow % v % n(c2)
        vel_max = max(vel_max, Flow % v % n(c2))
        vel_min = min(vel_min, Flow % v % n(c2))
      end if
    end if
  end do

  call Global % Sum_Real(area_in)
  call Global % Sum_Real(vol_in)
  vel_in = vol_in / area_in

  ! if(First_Proc()) then
  !   print *, '@ User_Mod_Beginning_Of_Time_Step'
  !   print *, '@ average inflow velocity = ', vel_in
  ! end if

  !------------------------------------------------!
  !   Set the average inflow velocity everywhere   !
  !------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if( Grid % Bnd_Cond_Type(c2) .eq. PRESSURE ) then
        Flow % v % n(c1) = vel_in
        Flow % v % n(c2) = vel_in
      end if
    end if
  end do

  end subroutine
