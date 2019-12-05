!==============================================================================!
  subroutine Swarm_Mod_Find_Nearest_Cell(swarm, k, n_parts_in_buffers)
!------------------------------------------------------------------------------!
!   Finds a cell closest to a particle.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k                   ! particle number
  integer                  :: n_parts_in_buffers
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: c, lc       ! cell, local cell
  integer                      :: cc, cb      ! closest cell and bnd. cell
  real                         :: xc, yc, zc  ! cell center coordinates
  real                         :: d_sq        ! distance squared
  real                         :: min_dc      ! minimum distance cell
  real                         :: min_db      ! minimum distance bnd cell
  real                         :: min_dc_glob
  real                         :: min_db_glob
!==============================================================================!

  ! Take aliases
  grid => swarm % pnt_grid
  part => swarm % particle(k)

  ! Initialize cc and cb to zero (important)
  cc = 0
  cb = 0

  !-----------------------------------------------------------!
  !                                                           !
  !   Closest node is known, browse through cells around it   !
  !     (You are entering this part during the simulation)    !
  !                                                           !
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
      d_sq = (xc - part % x_n)**2  &
           + (yc - part % y_n)**2  &
           + (zc - part % z_n)**2

      ! Finding the closest cell of those 8 neighbors
      if(c > 0) then
        if(d_sq < min_dc) then
          min_dc = d_sq  ! new minimum distance
          cc = c         ! cc is the closest inside cell index
        end if
      end if

      if(c < 0) then
        if(d_sq < min_db) then
          min_db = d_sq  ! new minimum distance
          cb = c         ! cb is the closest boundary cell index
        end if
      end if

    end do  ! browse through cells of closest node

    !--------------------------------------!
    !   Save the last value of cc and cb   !
    !     (the indices of closest cells)   !
    !--------------------------------------!
    part % cell     = cc
    part % bnd_cell = cb  ! this will be zero if no boundary cells
    part % bnd_face =  0

    ! If inside cell is closer than boundary ...
    ! ... cell don't store the boundary cell
    if(min_dc < min_db) then
      part % bnd_cell = 0
      part % bnd_face = 0
    end if
    if(min_db < min_dc) then
      if(cb .ne. 0) part % bnd_face = grid % cells_bnd_face(cb)
    end if

    !-----------------------------------------!
    !   Detect if particle entered a buffer   !
    !-----------------------------------------!
    part % buff = part % proc  ! assume buffer was not entered

    ! If processor number in the cell is differnt than this_proc 
    ! (and part % proc in this case) you entered the buffer
    if(grid % comm % cell_proc(cc) .ne. part % proc) then
      part % buff = grid % comm % cell_proc(cc)  ! store buffer process
      n_parts_in_buffers = n_parts_in_buffers + 1
    end if

  !---------------------------------------------------------!
  !                                                         !
  !   Closest node is not known, browse through all cells   !
  !             (You will be entering here only             !
  !             after the particle is introduced)           !
  !                                                         !
  !---------------------------------------------------------!
  else

    min_dc = HUGE  ! initialize minimum distance to cells inside
    min_db = HUGE  ! initialize minimum distance to boundary cells

    !-----------------------------!
    !   Cells inside the domain   !
    !-----------------------------!
    do c = 1, grid % n_cells - grid % comm % n_buff_cells

      ! Take cell centre
      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      ! Distance squared from the particle to cell centre
      d_sq = (xc - part % x_n)**2  &
           + (yc - part % y_n)**2  &
           + (zc - part % z_n)**2

      if(d_sq < min_dc) then
        min_dc = d_sq  ! new minimum distance
        cc = c         ! cc is the closest inside cell index
      end if

    end do  ! browswe through all cells

    part % cell = cc
    min_dc_glob = min_dc
    if(n_proc > 1) then
      call Comm_Mod_Global_Min_Real(min_dc_glob)
    end if

    part % proc = 0
    part % buff = 0
    if( min_dc .eq. min_dc_glob ) then
      part % proc = this_proc
      part % buff = this_proc
    end if

    !--------------------!
    !   Boundary cells   !
    !--------------------!
    do c = -grid % n_bnd_cells, -1

      ! Take cell centre
      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      ! Distance squared from the particle to cell centre
      d_sq = (xc - part % x_n)**2  &
           + (yc - part % y_n)**2  &
           + (zc - part % z_n)**2

      if(d_sq < min_db) then
        min_db = d_sq  ! new minimum distance
        cb = c         ! cc is the closest boundary cell index
      end if

    end do  ! browswe through all cells

    !--------------------------------------!
    !   Save the last value of cc and cb   !
    !     (the indices of closest cells)   !
    !--------------------------------------!
    part % bnd_cell = cb  ! this will be zero if no boundary cells here
    part % bnd_face =  0

    !--------------------------------------------!
    !   Check if particle is in this processor   !
    !--------------------------------------------!
    min_db_glob = min_db
    if(n_proc > 1) then
      call Comm_Mod_Global_Min_Real(min_db_glob)
    end if
    if(min_db .eq. min_db_glob .and. min_db < min_dc) then
      part % bnd_face = grid % cells_bnd_face(cb)
    end if

  end if  ! closest node is (not) known

  end subroutine
