!==============================================================================!
  subroutine Swarm_Mod_Print_Statistics(Swarm)
!------------------------------------------------------------------------------!
!   Prints particle statistics (still in early evolutionary stage)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: Part
  integer                      :: k, c, n_dep, n_ref, n_esc
  real                         :: avg_cfl, avg_re, avg_st
  real                         :: max_cfl, max_re, max_st
  character(DL)                :: line
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: T=38  ! indent
!==============================================================================!

  ! Take aliases for the Swarm
  Grid => Swarm % pnt_grid

  !----------------------------------------------!
  !   Sum and average data over all processors   !
  !----------------------------------------------!
  avg_cfl = 0
  avg_re  = 0
  avg_st  = 0
  max_cfl = -HUGE
  max_re  = -HUGE
  max_st  = -HUGE
  do k = 1, Swarm % n_particles
    Part => Swarm % Particle(k)
    if(Part % proc .eq. this_proc) then
      avg_cfl = avg_cfl + Part % cfl
      avg_re  = avg_re  + Part % re
      avg_st  = avg_st  + Part % st
      max_cfl = max(max_cfl, Part % cfl)
      max_re  = max(max_re,  Part % re)
      max_st  = max(max_st,  Part % st)
    end if
  end do
  call Comm_Mod_Global_Sum_Real(avg_cfl)
  call Comm_Mod_Global_Sum_Real(avg_st)
  call Comm_Mod_Global_Sum_Real(avg_re)
  call Comm_Mod_Global_Max_Real(max_cfl)
  call Comm_Mod_Global_Max_Real(max_re)
  call Comm_Mod_Global_Max_Real(max_st)
  avg_cfl = avg_cfl / real(Swarm % n_particles)
  avg_re  = avg_re  / real(Swarm % n_particles)
  avg_st  = avg_st  / real(Swarm % n_particles)

  n_dep = 0
  n_esc = 0
  n_ref = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % comm % cell_proc(c) .eq. this_proc) then  ! avoid buffer cells
      n_dep = n_dep + nint(Swarm % n_deposited(c))
      n_esc = n_esc + nint(Swarm % n_escaped(c))
      n_ref = n_ref + nint(Swarm % n_reflected(c))
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
    write(line(37+T:42+T),'(i6)') Swarm % n_particles
    print *, trim(line)
    line( 1+T:52+T) = ' #  Number of active particles    :               #'
    write(line(37+T:42+T),'(i6)') Swarm % n_particles - n_dep - n_esc
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
