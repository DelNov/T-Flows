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
  integer                      :: k, n_dep, n_ref, n_esc
  real                         :: avg_cfl, avg_re, avg_st
  real                         :: max_cfl, max_re, max_st
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !----------------------------------------------!
  !   Sum and average data over all processors   !
  !----------------------------------------------!
  avg_cfl = 0
  avg_re  = 0
  avg_st  = 0
  max_cfl = -HUGE
  max_re  = -HUGE
  max_st  = -HUGE
  do k = 1, swarm % n_particles
    part => swarm % particle(k)
    if(part % proc .eq. this_proc) then
      avg_cfl = avg_cfl + part % cfl
      avg_re  = avg_re  + part % re
      avg_st  = avg_st  + part % st
      max_cfl = max(max_cfl, part % cfl)
      max_re  = max(max_re,  part % re)
      max_st  = max(max_st,  part % st)
    end if
  end do
  call Comm_Mod_Global_Sum_Real(avg_cfl)
  call Comm_Mod_Global_Sum_Real(avg_st)
  call Comm_Mod_Global_Sum_Real(avg_re)
  call Comm_Mod_Global_Max_Real(max_cfl)
  call Comm_Mod_Global_Max_Real(max_re)
  call Comm_Mod_Global_Max_Real(max_st)
  avg_cfl = avg_cfl / swarm % n_particles
  avg_re  = avg_re  / swarm % n_particles
  avg_st  = avg_st  / swarm % n_particles

  n_dep = nint(sum(swarm % n_deposited(:)))
  n_esc = nint(sum(swarm % n_escaped(:)))
  n_ref = nint(sum(swarm % n_reflected(:)))
  call Comm_Mod_Global_Sum_Int(n_dep)
  call Comm_Mod_Global_Sum_Int(n_esc)
  call Comm_Mod_Global_Sum_Int(n_ref)

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  if(this_proc < 2) then
    write(*,'(a)') ' #================================================#'
    write(*,'(a)') ' #                   Swarm statistics'
    write(*,'(a)') ' #------------------------------------------------#'
    write(*,'(a,i6)') ' #  Total number of particles     : ',  &
                      swarm % n_particles
    write(*,'(a,i6)') ' #  Number of active particles    : ',  &
                      swarm % n_particles - n_dep - n_esc
    write(*,'(a,i6)')      ' #  Number of deposited particles : ',  n_dep
    write(*,'(a,i6)')      ' #  Number of escaped particles   : ',  n_esc
    write(*,'(a,1pe12.1)') ' #  Total number of reflections   : ',  real(n_ref)
    write(*,'(a)') ' #================================================#'
    write(*,'(a)') ' #     Characteristic non-dimensional numbers     #'
    write(*,'(a)') ' #          (average and maximum values)          #'
    write(*,'(a)') ' #------------------------------------------------#'
    write(*,'(a,2(1pe11.3))') ' #    Courant number   :',  avg_cfl, max_cfl
    write(*,'(a,2(1pe11.3))') ' #    Reynolds number  :',  avg_re,  max_re
    write(*,'(a,2(1pe11.3))') ' #    Stokes number    :',  avg_st,  max_st
    write(*,'(a)') ' #------------------------------------------------#'
  end if

  end subroutine
