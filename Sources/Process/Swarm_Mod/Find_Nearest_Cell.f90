!==============================================================================!
  subroutine Swarm_Mod_Find_Nearest_Cell(swarm, k)
!------------------------------------------------------------------------------!
!   Finds a cell closest to a particle.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: c, lc       ! cell, local cell
  integer                      :: cc, cb      ! closest cell and bnd. cell
  real                         :: xc, yc, zc  ! cell center coordinates
  real                         :: dc_sq       ! distance squared
  real                         :: min_dc      ! minimum distance cell
  real                         :: min_db      ! minimum distance bnd cell
!==============================================================================!

  ! Take aliases
  grid => swarm % pnt_grid
  part => swarm % particle(k)

  cc = 0
  cb = 0

  !-----------------------------------------------------------!
  !   Closest node is known, browse through cells around it   !
  !-----------------------------------------------------------!
  if(part % node .ne. 0) then

    min_dc = HUGE  ! initialize minimum distance to cells inside
    min_db = HUGE  ! initialize minimum distance to boundary cells

    do lc = 1, grid % nodes_n_cells(part % node)   ! local cell number
      c = grid % nodes_c(lc, part % node)          ! global cell number

      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      ! Distance squared from the particle to cell centre
      dc_sq = (xc - part % x)**2  &
            + (yc - part % y)**2  &
            + (zc - part % z)**2

      ! Finding the closest cell of those 8 neighbors
      if(c > 0) then
        if(dc_sq < min_dc) then
          min_dc = dc_sq  ! new minimum distance
          cc = c          ! cc is the closest inside cell index
        end if
      end if

      if(c < 0) then
        if(dc_sq < min_db) then
          min_db = dc_sq  ! new minimum distance
          cb = c          ! cb is the closest boundary cell index
        end if
      end if

    end do

  !---------------------------------------------------------!
  !   Closest node is not known, browse through all cells   !
  !---------------------------------------------------------!
  else

    min_dc = HUGE  ! initialize minimum distance to cells inside
    min_db = HUGE  ! initialize minimum distance to boundary cells

    do c = -grid % n_bnd_cells, grid % n_cells

      if(c .ne. 0) then

        ! Take cell centre
        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

        ! Distance squared from the particle to cell centre
        dc_sq = (xc - part % x)**2  &
              + (yc - part % y)**2  &
              + (zc - part % z)**2

        if(c > 0) then
          if(dc_sq < min_dc) then
            min_dc = dc_sq  ! new minimum distance
            cc = c          ! cc is the closest inside cell index
          end if
        end if

        if(c < 0) then
          if(dc_sq < min_db) then
            min_db = dc_sq  ! new minimum distance
            cb = c          ! cc is the closest boundary cell index
          end if
        end if

      end if  ! c .ne. 0

    end do

  end if  ! closest node is (not) known

  !--------------------------------------!
  !   Save the last value of cc and cb   !
  !     (the indices of closest cells)   !
  !--------------------------------------!
  part % cell     = cc
  part % bnd_cell = cb
  part % bnd_face =  0
  if(cb < 0) part % bnd_face = grid % cells_bnd_face(cb)

  ! If inside cell is closer than boundary ...
  ! ... cell don't store the boundary cell
  if(min_dc < min_db) then
    part % bnd_cell = 0
    part % bnd_face = 0
  end if

  end subroutine
