!==============================================================================!
  subroutine Find_Nearest_Cell(Point, n_parts_in_buffers, locally)
!------------------------------------------------------------------------------!
!>  This subroutine determines the closest cell to a given point within the
!>  computational grid. It is a critical function for accurately positioning
!>  points, such as particles, within the grid space, and for managing their
!>  interactions with grid elements during simulations.  It also estimates if
!>  the point traverses inter-processor boundaries, entering into buffers.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Prepares necessary variables for distance computation    !
!     and identifies whether the search is local to a subdomain or global.     !
!   * Proximity analysis: Determines the closest cell to the point by          !
!     calculating the distance to cell centers, both for internal and          !
!     boundary cells.                                                          !
!   * Closest cell identification: Identifies the nearest cell and boundary    !
!     cell based on the computed distances.                                    !
!   * Buffer zone detection: Establishes if the point has entered a buffer     !
!     zone, crucial for parallel processing scenarios.                         !
!   * Global minimum distance: In parallel runs, finds the global minimum      !
!     distance to ensure the point is in the correct processor's subdomain.    !
!   * Point update: Updates the point's properties with the identified closest !
!     cell and processor information.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type), target :: Point               !! point object
  integer                   :: n_parts_in_buffers  !! counter for the number of
                                                   !! particles in buffer zones
  logical, optional         :: locally             !! flag to restrict search
                                                   !! to the local subdomain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, i_cel    ! cell, local cell
  integer                  :: cc, cb      ! closest cell and bnd. cell
  logical                  :: local
  real                     :: xc, yc, zc  ! cell center coordinates
  real                     :: d_sq        ! distance squared
  real                     :: min_dc      ! minimum distance cell
  real                     :: min_db      ! minimum distance bnd cell
  real                     :: min_dc_glob
  real                     :: min_db_glob
!==============================================================================!

  ! Take aliases
  Grid => Point % pnt_grid

  local = .false.
  if(present(locally)) local = locally

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

    ! If processor number in the cell is differnt than This_Proc()
    ! (and Point % proc in this case) you entered the buffer
    if(Grid % Comm % cell_proc(cc) .ne. Point % proc) then
      Point % buff = Grid % Comm % cell_proc(cc)  ! store buffer process
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
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells

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
    if(Parallel_Run() .and. .not. local) then
      call Global % Min_Real(min_dc_glob)
    end if

    Point % proc = 0
    Point % buff = 0
    if( min_dc .eq. min_dc_glob ) then
      Point % proc = This_Proc()
      Point % buff = This_Proc()
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
    if(Parallel_Run() .and. .not. local) then
      call Global % Min_Real(min_db_glob)
    end if

  end if  ! closest node is (not) known

  !------------------------------!
  !                              !
  !   Check if the particle is   !
  !    really inside the cell    !
  !     (Didn't really work)     !
  !                              !
  !------------------------------!
  ! if(Point % proc .eq. This_Proc()                      &
  !    .and. .not. Grid % Is_Point_In_Cell(Point % cell,  &
  !                                        Point % x_n,   &
  !                                        Point % y_n,   &
  !                                        Point % z_n)) then
  !   ! Mark it as -1, like invalid
  !   Point % cell     = -1
  !   Point % bnd_cell =  0
  ! end if

  end subroutine
