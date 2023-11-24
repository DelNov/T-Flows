!==============================================================================!
  subroutine Check_Periodicity(Swarm, k, n_parts_in_buffers)
!------------------------------------------------------------------------------!
!> This subroutine checks if a particle in a swarm has crossed a periodic
!> boundary in the computational domain. It is crucial for maintaining
!> accurate tracking of particles in simulations with periodic boundary
!> conditions, ensuring particles are correctly relocated to their new
!> positions within the domain.
!------------------------------------------------------------------------------!
! Functionality                                                                !
!                                                                              !
! * Periodic boundary handling: Identifies if particles have crossed periodic  !
!   boundaries and adjusts their positions accordingly.                        !
! * Particle relocation: Ensures particles that cross a periodic boundary are  !
!   accurately relocated to their corresponding shadow positions.              !
! * Buffer zone management: Marks particles that move into buffer zones,       !
!   aiding in parallel processing scenarios.                                   !
! * Distance calculation: Computes distances to ascertain the closest shadow   !
!   position of a particle after crossing a boundary.                          !
! * Cell update: Updates the cell and node information of a particle after     !
!   crossing a periodic boundary.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm               !! the swarm of particles
  integer                   :: k                   !! particle number (rank)
  integer                   :: n_parts_in_buffers  !! counter for the number of
                                                   !! particles in buffer zones
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: Part
  integer                      :: s, c1, c2
  real                         :: xc1, xc2, yc1, yc2, zc1, zc2, d_sq_1, d_sq_2
!==============================================================================!

  ! Take aliases
  Grid => Swarm % pnt_grid

  !-------------------------------------!
  !                                     !
  !   Mark cells which have particles   !
  !                                     !
  !-------------------------------------!
  Swarm % cell_has_particles(:) = .false.  ! assume no cells has particles

  Part => Swarm % Particle(k)
  if(Part % proc .eq. This_Proc()) then
    Swarm % cell_has_particles(Part % cell) = .true.
  end if

  !-----------------------------------!
  !                                   !
  !   Browse through periodic faces   !
  !                                   !
  !-----------------------------------!
  do s = 1, Grid % n_faces

    if(Grid % faces_s(s) .ne. 0) then

      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      !------------------------------!
      !   If c1 contains particles   !
      !------------------------------!
      if( Swarm % cell_has_particles(c1) ) then

        Part => Swarm % Particle(k)

        if(Part % proc .eq. This_Proc()) then

          xc1 = Grid % xc(c1)
          yc1 = Grid % yc(c1)
          zc1 = Grid % zc(c1)

          ! Take position of c2 at its shadow position
          xc2 = Grid % xc(c1) + Grid % dx(s)
          yc2 = Grid % yc(c1) + Grid % dy(s)
          zc2 = Grid % zc(c1) + Grid % dz(s)

          d_sq_1 = (xc1 - Part % x_n)**2  &
                 + (yc1 - Part % y_n)**2  &
                 + (zc1 - Part % z_n)**2

          d_sq_2 = (xc2 - Part % x_n)**2  &
                 + (yc2 - Part % y_n)**2  &
                 + (zc2 - Part % z_n)**2

          ! Particle is indeed closer to c2's shadow, mark it as "lost"
          if( d_sq_2 < d_sq_1 ) then

            ! Move c1 to its shadow position
            xc1 = Grid % xc(c2) - Grid % dx(s)
            yc1 = Grid % yc(c2) - Grid % dy(s)
            zc1 = Grid % zc(c2) - Grid % dz(s)

            Part % x_n = xc1 + Part % x_n - Grid % xc(c1)
            Part % y_n = yc1 + Part % y_n - Grid % yc(c1)
            Part % z_n = zc1 + Part % z_n - Grid % zc(c1)
            Part % x_o = xc1 + Part % x_o - Grid % xc(c1)
            Part % y_o = yc1 + Part % y_o - Grid % yc(c1)
            Part % z_o = zc1 + Part % z_o - Grid % zc(c1)
                              !<---- dist(Part, c1) ---->!
            Part % cell = c2
            ! call Swarm_Mod_Find_Nearest_Node(Swarm, k)
            call Part % Find_Nearest_Node()

            ! If c2 is in the buffer, tell that particle wants to go there
            if(Grid % Comm % cell_proc(c2) .ne.  &
               Grid % Comm % cell_proc(c1)) then
              Part % buff = Grid % Comm % cell_proc(c2)
              n_parts_in_buffers = n_parts_in_buffers + 1
            end if
          end if
        end if

      end if    ! c1 has particles

      !------------------------------!
      !   If c2 contains particles   !
      !------------------------------!
      if( Swarm % cell_has_particles(c2) ) then

        Part => Swarm % Particle(k)

        if(Part % proc .eq. This_Proc()) then

          ! Take position of c1 at its shadow position
          xc1 = Grid % xc(c2) - Grid % dx(s)
          yc1 = Grid % yc(c2) - Grid % dy(s)
          zc1 = Grid % zc(c2) - Grid % dz(s)

          xc2 = Grid % xc(c2)
          yc2 = Grid % yc(c2)
          zc2 = Grid % zc(c2)

          d_sq_1 = (xc1 - Part % x_n)**2  &
                 + (yc1 - Part % y_n)**2  &
                 + (zc1 - Part % z_n)**2

          d_sq_2 = (xc2 - Part % x_n)**2  &
                 + (yc2 - Part % y_n)**2  &
                 + (zc2 - Part % z_n)**2

          ! Particle is indeed closer to c1's shadow, mark it as "lost"
          if( d_sq_1 < d_sq_2 ) then

            ! Move c2 to its shadow position
            xc2 = Grid % xc(c1) + Grid % dx(s)
            yc2 = Grid % yc(c1) + Grid % dy(s)
            zc2 = Grid % zc(c1) + Grid % dz(s)

            Part % x_n = xc2 + Part % x_n - Grid % xc(c2)
            Part % y_n = yc2 + Part % y_n - Grid % yc(c2)
            Part % z_n = zc2 + Part % z_n - Grid % zc(c2)
            Part % x_o = xc2 + Part % x_o - Grid % xc(c2)
            Part % y_o = yc2 + Part % y_o - Grid % yc(c2)
            Part % z_o = zc2 + Part % z_o - Grid % zc(c2)
                              !<---- dist(Part, c2) ---->!
            Part % cell = c1
            ! call Swarm_Mod_Find_Nearest_Node(Swarm, k)
            call Part % Find_Nearest_Node()

            ! If c1 is in the buffer, tell that particle wants to go there
            if(Grid % Comm % cell_proc(c1) .ne.  &
               Grid % Comm % cell_proc(c2)) then
              Part % buff = Grid % Comm % cell_proc(c1)
              n_parts_in_buffers = n_parts_in_buffers + 1
            end if
          end if
        end if

      end if    ! c2 has particles

    end if  ! is a shadow face

  end do

  end subroutine
