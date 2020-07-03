!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer, intent(in)           :: n         ! time step
  integer, intent(in)           :: n_stat_t  ! 1st step for turb. stat.
  integer, intent(in)           :: n_stat_p  ! 1st step for swarm stat.
  real,    intent(in)           :: time      ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: k, n_parts_in_buffers
  real    :: dx
!==============================================================================!

  !----------------------!
  !   2nd time step on   !
  !----------------------!
  if(n .gt.  15) then     ! should be after the flow is developed
    call Swarm_Mod_Advance_Particles(swarm, n, n_stat_p)
  end if

  if(this_proc < 2) then
    write(*,'(a,i4,a,i4,a,i4)')                        &
             "# trapped particles: ",  swarm % cnt_d,  &
             " escaped particles: ",   swarm % cnt_e,  &
             " reflected particles: ", swarm % cnt_r
  end if

  end subroutine
