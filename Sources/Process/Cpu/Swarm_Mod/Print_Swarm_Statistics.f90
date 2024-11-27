!==============================================================================!
  subroutine Print_Swarm_Statistics(Swarm)
!------------------------------------------------------------------------------!
!>  Gathers and prints statistical data on the particles within a swarm. This
!>  subroutine is an integral part of monitoring and analyzing the behavior of
!>  particles in computational fluid dynamics simulations, providing insights
!>  into their collective dynamics and individual properties.
!------------------------------------------------------------------------------!
!   Outputs and Calculations                                                   !
!                                                                              !
!   * Particle counts: Reports total, active, deposited, and escaped           !
!     particles, offering a snapshot of particle distribution in the domain.   !
!   * Reflection data: Quantifies the number of particle reflections,          !
!     indicating interactions with boundaries.                                 !
!   * Non-dimensional numbers: Presents average and maximum values for         !
!     Courant, Reynolds, and Stokes numbers, crucial for assessing particle    !
!     dynamics relative to the fluid flow.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: Part
  integer                      :: k, c, n_dep, n_ref, n_esc
  real                         :: avg_cfl, avg_re, avg_st
  real                         :: max_cfl, max_re, max_st
  character(DL)                :: line
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: T=39  ! indent
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
    if(Part % proc .eq. This_Proc()) then
      avg_cfl = avg_cfl + Part % cfl
      avg_re  = avg_re  + Part % re
      avg_st  = avg_st  + Part % st
      max_cfl = max(max_cfl, Part % cfl)
      max_re  = max(max_re,  Part % re)
      max_st  = max(max_st,  Part % st)
    end if
  end do
  call Global % Sum_Real(avg_cfl)
  call Global % Sum_Real(avg_st)
  call Global % Sum_Real(avg_re)
  call Global % Max_Real(max_cfl)
  call Global % Max_Real(max_re)
  call Global % Max_Real(max_st)
  avg_cfl = avg_cfl / real(Swarm % n_particles)
  avg_re  = avg_re  / real(Swarm % n_particles)
  avg_st  = avg_st  / real(Swarm % n_particles)

  n_dep = 0
  n_esc = 0
  n_ref = 0
  do c = -Grid % n_bnd_cells, -1
    if(Cell_In_This_Proc(c)) then  ! avoid buffer cells
      n_dep = n_dep + nint(Swarm % n_deposited(c))
      n_esc = n_esc + nint(Swarm % n_escaped(c))
      n_ref = n_ref + nint(Swarm % n_reflected(c))
    end if
  end do
  call Global % Sum_Int(n_dep)
  call Global % Sum_Int(n_esc)
  call Global % Sum_Int(n_ref)

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  if(First_Proc()) then
    line( 1:160) = ' '
    line( 1+T:52+T) = ' #================================================#'
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #                Swarm statistics                #'
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #------------------------------------------------#'
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #  Total number of particles     :               #'
    write(line(37+T:42+T),'(i6)') Swarm % n_particles
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #  Number of active particles    :               #'
    write(line(37+T:42+T),'(i6)') Swarm % n_particles - n_dep - n_esc
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #  Number of deposited particles :               #'
    write(line(37+T:42+T),'(i6)') n_dep
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #  Number of escaped particles   :               #'
    write(line(37+T:42+T),'(i6)') n_esc
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #  Total number of reflections   :               #'
    write(line(37+T:48+T),'(1pe12.1)') real(n_ref)
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #================================================#'
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #     Characteristic non-dimensional numbers     #'
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #          (average and maximum values)          #'
    print '(a)', trim(line)
    line( 1+T:52+T) = ' #------------------------------------------------#'
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #    Courant number   :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_cfl, max_cfl
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #    Reynolds number  :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_re,  max_re
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #    Stokes number    :                          #'
    write(line(26+T:48+T),'(2(1pe11.3))') avg_st,  max_st
    print '(a)', trim(line)

    line( 1+T:52+T) = ' #------------------------------------------------#'
    print '(a)', trim(line)
  end if

  end subroutine
