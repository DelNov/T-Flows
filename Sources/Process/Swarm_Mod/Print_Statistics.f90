!==============================================================================!
  subroutine Swarm_Mod_Print_Statistics(swarm)
!------------------------------------------------------------------------------!
!   Prints particle statistics (still in early evolutionary stage)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: k
  real                         :: avg_part_cfl, avg_part_re, avg_part_st
  real                         :: max_part_cfl, max_part_re, max_part_st
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  avg_part_cfl = 0
  avg_part_re  = 0
  avg_part_st  = 0
  max_part_cfl = -HUGE
  max_part_re  = -HUGE
  max_part_st  = -HUGE
  do k = 1, swarm % n_particles
    part => swarm % particle(k)
    if(part % proc .eq. this_proc) then
      avg_part_cfl = avg_part_cfl + part % cfl
      avg_part_re  = avg_part_re  + part % re
      avg_part_st  = avg_part_st  + part % st
      max_part_cfl = max(max_part_cfl, part % cfl)
      max_part_re  = max(max_part_re,  part % re)
      max_part_st  = max(max_part_st,  part % st)
    end if
  end do
  call Comm_Mod_Global_Sum_Real(avg_part_cfl)
  call Comm_Mod_Global_Sum_Real(avg_part_st)
  call Comm_Mod_Global_Sum_Real(avg_part_re)
  call Comm_Mod_Global_Max_Real(max_part_cfl)
  call Comm_Mod_Global_Max_Real(max_part_re)
  call Comm_Mod_Global_Max_Real(max_part_st)
  avg_part_cfl = avg_part_cfl / swarm % n_particles
  avg_part_re  = avg_part_re  / swarm % n_particles
  avg_part_st  = avg_part_st  / swarm % n_particles
  if(this_proc < 2) then
    write(*, '(a)') ' #====================================================='
    write(*, '(a)') ' #                   Swarm statistics'
    write(*, '(a)') ' #-----------------------------------------------------'
    write(*, '(a,i7)') ' # Total number of particles     : ',  &
                       swarm % n_particles
    write(*, '(a,i7)') ' # Number of active particles    : ',  &
                       swarm % n_particles - nint(sum(swarm % n_deposited(:)))
    write(*, '(a,i7)') ' # Number of deposited particles : ',  &
                       nint(sum(swarm % n_deposited(:)))
    write(*, '(a,i7)') ' # Number of escaped particles   : ',  &
                       nint(sum(swarm % n_escaped(:)))
    write(*, '(a,1pe13.1)') ' # Total number of reflections   : ',  &
                       sum(swarm % n_reflected(:))
    write(*, '(a)') ' #-----------------------------------------------------'
    write(*,'(a,2(1pe11.3))') ' # Average and maximum Courant number  : ',  &
                          avg_part_cfl, max_part_cfl
    write(*,'(a,2(1pe11.3))') ' # Average and maximum Reynolds number : ',  &
                          avg_part_re,  max_part_re
    write(*,'(a,2(1pe11.3))') ' # Average and maximum Stokes number   : ',  &
                          avg_part_st,  max_part_st
    write(*, '(a)') ' #-----------------------------------------------------'
  end if

  end subroutine
