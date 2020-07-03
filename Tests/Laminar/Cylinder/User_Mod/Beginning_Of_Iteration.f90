!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, vof
  integer                  :: s, c1, c2
  real                     :: vol_in, area_in, vel_in
  real                     :: vel_max, vel_min
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  !-------------------------------------!
  !   Compute average inflow velocity   !
  !-------------------------------------!
  vol_in  = 0.0
  area_in = 0.0

  vel_max = -HUGE
  vel_min = +HUGE
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE ) then
        area_in = area_in + grid % s(s)
        vol_in  = vol_in  + grid % s(s) * flow % v % n(c2)
        vel_max = max(vel_max, flow % v % n(c2))
        vel_min = min(vel_min, flow % v % n(c2))
      end if
    end if
  end do

  call Comm_Mod_Global_Sum_Real(area_in)
  call Comm_Mod_Global_Sum_Real(vol_in)
  vel_in = vol_in / area_in

  ! if(this_proc < 2) then
  !   print *, '@ User_Mod_Beginning_Of_Time_Step'
  !   print *, '@ average inflow velocity = ', vel_in
  ! end if

  !------------------------------------------------!
  !   Set the average inflow velocity everywhere   !
  !------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE ) then
        flow % v % n(c1) = vel_in
        flow % v % n(c2) = vel_in
      end if
    end if
  end do

  end subroutine
