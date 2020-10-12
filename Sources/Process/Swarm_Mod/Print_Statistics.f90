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
  integer                      :: k, c, n_dep, n_ref, n_esc
  real                         :: avg_cfl, avg_re, avg_st
  real                         :: max_cfl, max_re, max_st
  character(DL)                :: line
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: T=38  ! indent
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

  n_dep = 0
  n_esc = 0
  n_ref = 0
  do c = -grid % n_bnd_cells, -1
    if(grid % comm % cell_proc(c) .eq. this_proc) then  ! avoid buffer cells
      n_dep = n_dep + nint(swarm % n_deposited(c))
      n_esc = n_esc + nint(swarm % n_escaped(c))
      n_ref = n_ref + nint(swarm % n_reflected(c))
    end if
  end do
  call Comm_Mod_Global_Sum_Int(n_dep)
  call Comm_Mod_Global_Sum_Int(n_esc)
  call Comm_Mod_Global_Sum_Int(n_ref)

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  if(this_proc < 2) then
    line( 1:160) = ' '
    line( 1+T:52+T) = ' #================================================#'
    print *, trim(line)
    line( 1+T:52+T) = ' #                   Swarm statistics             #'
    print *, trim(line)
    line( 1+T:52+T) = ' #------------------------------------------------#'
    print *, trim(line)

    line( 1+T:52+T) = ' #  Total number of particles     :               #'
    write(line(37+T:42+T),'(i6)') swarm % n_particles
    print *, trim(line)
    line( 1+T:52+T) = ' #  Number of active particles    :               #'
    write(line(37+T:42+T),'(i6)') swarm % n_particles - n_dep - n_esc
    print *, trim(line)

    line( 1+T:52+T) = ' #  Number of deposited particles :               #'
    write(line(37+T:42+T),'(i6)') n_dep
    print *, trim(line)
    line( 1+T:52+T) = ' #  Number of escaped particles   :               #'
    write(line(37+T:42+T),'(i6)') n_esc
    print *, trim(line)
    line( 1+T:52+T) = ' #  Total number of reflections   :               #'
    write(line(37+T:48+T),'(1pe12.1)') real(n_ref)
    print *, trim(line)

    line( 1+T:52+T) = ' #================================================#'
    print *, trim(line)
    line( 1+T:52+T) = ' #     Characteristic non-dimensional numbers     #'
    print *, trim(line)
    line( 1+T:52+T) = ' #          (average and maximum values)          #'
    print *, trim(line)
    line( 1+T:52+T) = ' #------------------------------------------------#'
    print *, trim(line)

    line( 1+T:52+T) = ' #    Courant number   :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_cfl, max_cfl
    print *, trim(line)

    line( 1+T:52+T) = ' #    Reynolds number  :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_re,  max_re
    print *, trim(line)

    line( 1+T:52+T) = ' #    Stokes number    :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_st,  max_st
    print *, trim(line)

    line( 1+T:52+T) = ' #------------------------------------------------#'
    print *, trim(line)
  end if

  end subroutine
