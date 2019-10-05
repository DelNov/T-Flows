!==============================================================================!
  subroutine Swarm_Mod_Check_Periodicity(swarm)
!------------------------------------------------------------------------------!
!   Check if particle left the cell (which means periodic domain as well)      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: k, ps, s, c1, c2
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid

  !-------------------------------------!
  !   Mark cells which have particles   !
  !-------------------------------------!
  swarm % cell_has_particles(:) = .false.  ! assume no cells has particles

  do k = 1, swarm % n_particles
    part => swarm % particle(k)
    if(part % proc .eq. this_proc) then
      swarm % cell_has_particles(part % cell) = .true.
    end if
  end do

  !-----------------------------------!
  !   Browse through periodic faces   !
  !-----------------------------------!
  do ps = 1, grid % n_per_faces
    s = grid % per_faces(ps)      ! take periodic face number

    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    ! If any of the surrounding cells holds particles, do the check
    if( swarm % cell_has_particles(c1) .or.   &
        swarm % cell_has_particles(c2) ) then

      ! To be done if particle crosses face s ...

    end if
  end do

  end subroutine
