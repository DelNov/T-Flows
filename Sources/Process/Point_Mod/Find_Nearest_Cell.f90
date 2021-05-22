!==============================================================================!
  subroutine Find_Nearest_Cell(Point, n_parts_in_buffers)
!------------------------------------------------------------------------------!
!   Finds a cell closest to the point                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type), target :: Point
  integer                   :: n_parts_in_buffers
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, i_cel    ! cell, local cell
  integer                  :: cc, cb      ! closest cell and bnd. cell
  real                     :: xc, yc, zc  ! cell center coordinates
  real                     :: d_sq        ! distance squared
  real                     :: min_dc      ! minimum distance cell
  real                     :: min_db      ! minimum distance bnd cell
  real                     :: min_dc_glob
  real                     :: min_db_glob
!==============================================================================!

  ! Take aliases
  Grid => Point % pnt_grid

  ! Initialize cc and cb to zero (important)
  cc = 0
  cb = 0

  !-----------------------------------------------------------!
  !                                                           !
  !   Closest node is known, browse through cells around it   !
  !     (You are entering this part during the simulation)    !
  !                                                           !
  !-----------------------------------------------------------!
  if(Point % node .ne. 0) then

    min_dc = HUGE  ! initialize minimum distance to cells inside
    min_db = HUGE  ! initialize minimum distance to boundary cells

    do i_cel = 1, Grid % nodes_n_cells(Point % node)   ! local cell number
      c = Grid % nodes_c(i_cel, Point % node)          ! global cell number

      xc = Grid % xc(c)
      yc = Grid % yc(c)
      zc = Grid % zc(c)

      ! Distance squared from the point to cell centre
      d_sq = (xc - Point % x_n)**2  &
           + (yc - Point % y_n)**2  &
           + (zc - Point % z_n)**2

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
    Point % cell     = cc
    Point % bnd_cell = cb  ! this will be zero if no boundary cells

    ! If inside cell is closer than boundary ...
    ! ... cell don't store the boundary cell
    if(min_dc < min_db) then
      Point % bnd_cell = 0
    end if

    !--------------------------------------!
    !   Detect if point entered a buffer   !
    !--------------------------------------!
    Point % buff = Point % proc  ! assume buffer was not entered

    ! If processor number in the cell is differnt than this_proc 
    ! (and Point % proc in this case) you entered the buffer
    if(Grid % comm % cell_proc(cc) .ne. Point % proc) then
      Point % buff = Grid % comm % cell_proc(cc)  ! store buffer process
      n_parts_in_buffers = n_parts_in_buffers + 1
    end if

  !---------------------------------------------------------!
  !                                                         !
  !   Closest node is not known, browse through all cells   !
  !             (You will be entering here only             !
  !              after the point is introduced)             !
  !                                                         !
  !---------------------------------------------------------!
  else

    min_dc = HUGE  ! initialize minimum distance to cells inside
    min_db = HUGE  ! initialize minimum distance to boundary cells

    !-----------------------------!
    !   Cells inside the domain   !
    !-----------------------------!
    do c = 1, Grid % n_cells - Grid % comm % n_buff_cells

      ! Take cell centre
      xc = Grid % xc(c)
      yc = Grid % yc(c)
      zc = Grid % zc(c)

      ! Distance squared from the point to cell centre
      d_sq = (xc - Point % x_n)**2  &
           + (yc - Point % y_n)**2  &
           + (zc - Point % z_n)**2

      if(d_sq < min_dc) then
        min_dc = d_sq  ! new minimum distance
        cc = c         ! cc is the closest inside cell index
      end if

    end do  ! browse through all cells

    Point % cell = cc
    min_dc_glob = min_dc
    if(n_proc > 1) then
      call Comm_Mod_Global_Min_Real(min_dc_glob)
    end if

    Point % proc = 0
    Point % buff = 0
    if( min_dc .eq. min_dc_glob ) then
      Point % proc = this_proc
      Point % buff = this_proc
    end if

    !--------------------!
    !   Boundary cells   !
    !--------------------!
    do c = -Grid % n_bnd_cells, -1

      ! Take cell centre
      xc = Grid % xc(c)
      yc = Grid % yc(c)
      zc = Grid % zc(c)

      ! Distance squared from the point to cell centre
      d_sq = (xc - Point % x_n)**2  &
           + (yc - Point % y_n)**2  &
           + (zc - Point % z_n)**2

      if(d_sq < min_db) then
        min_db = d_sq  ! new minimum distance
        cb = c         ! cc is the closest boundary cell index
      end if

    end do  ! browswe through all cells

    !--------------------------------------!
    !   Save the last value of cc and cb   !
    !     (the indices of closest cells)   !
    !--------------------------------------!
    Point % bnd_cell = cb  ! this will be zero if no boundary cells here

    !-----------------------------------------!
    !   Check if point is in this processor   !
    !-----------------------------------------!
    min_db_glob = min_db
    if(n_proc > 1) then
      call Comm_Mod_Global_Min_Real(min_db_glob)
    end if

  end if  ! closest node is (not) known

  end subroutine
