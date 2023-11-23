!==============================================================================!
  subroutine Distribute_Ranges(Dom, Grid)
!------------------------------------------------------------------------------!
!>  Distribute ranges (defined in .dom file) to specific boundary conditions
!>  in a computational grid.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Allocate regions in grid: The subroutine starts by allocating memory     !
!     for regions within the Grid object, based on the number of ranges in Dom !
!   * Initialization: Initializes the number of boundary conditions to zero.   !
!   * Iterating over ranges: Loops through each range in Dom, performing the   !
!     following steps:                                                         !
!     - Retrieves the associated block and its resolution.                     !
!     - Determines the default values and boundary faces based on the          !
!       range's specifications.                                                !
!     - Differentiates between boundary conditions prescribed with mnemonics   !
!       (like 'IMIN', 'IMAX') and those prescribed explicitly by coordinates.  !
!   * Storing Boundary Conditions:                                             !
!     - Checks if a boundary condition with the same name already exists       !
!     - If not, adds the new boundary condition to Grid.                       !
!     - Assigns a negative index to the corresponding cells in the grid,       !
!       indicating the boundary condition applied to those cells.              !
!   * Updating Grid Information: Finally, stores the total number of           !
!     boundary conditions in Grid and prints the regions list.                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom   !! domain in which the grid is being generated
  type(Grid_Type)     :: Grid  !! grid being generated
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, i, j, k, n, c, r
  integer :: n_bnd                         ! number of boundary conditions
  integer :: is, js, ks, ie, je, ke, face
  integer :: ci, cj, ck
  logical :: found
!==============================================================================!

  !-----------------------------------------!
  !   Insertion of the boundary condition   !
  !        and materials information        !
  !-----------------------------------------!

  ! This is too much memory but that's OK
  !  (+1 is to store the default values)
  call Grid % Allocate_Regions(Dom % n_ranges)

  ! Initialize number of boundary conditions
  n_bnd = 0

  do n = 1, Dom % n_ranges

    b = Dom % ranges(n) % block

    ! Block resolution
    ci = Dom % blocks(b) % resolutions(1)-1
    cj = Dom % blocks(b) % resolutions(2)-1
    ck = Dom % blocks(b) % resolutions(3)-1

    ! Default values
    is = 1
    ie = ci
    js = 1
    je = cj
    ks = 1
    ke = ck

    ! Boundary conditions prescribed with mnemonics
    if(Dom % ranges(n) % face .eq. 'IMIN') then
      ie   = 1
      face = 5
    else if(Dom % ranges(n) % face .eq. 'IMAX') then
      is   = ci
      face = 3
    else if(Dom % ranges(n) % face .eq. 'JMIN') then
      je   = 1
      face = 2
    else if(Dom % ranges(n) % face .eq. 'JMAX') then
      js   = cj
      face = 4
    else if(Dom % ranges(n) % face .eq. 'KMIN') then
      ke   = 1
      face = 1
    else if(Dom % ranges(n) % face .eq. 'KMAX') then
      ks   = ck
      face = 6

    ! Boundary conditions  prescribed explicitly
    !  (error prone and difficult, but might be usefull)
    else
      is = Dom % ranges(n) % is
      js = Dom % ranges(n) % js
      ks = Dom % ranges(n) % ks
      ie = Dom % ranges(n) % ie
      je = Dom % ranges(n) % je
      ke = Dom % ranges(n) % ke
      face = 0
      if( (is .eq. ie).and.(is .eq.  1) ) face = 5
      if( (is .eq. ie).and.(is .eq. ci) ) face = 3
      if( (js .eq. je).and.(js .eq.  1) ) face = 2
      if( (js .eq. je).and.(js .eq. cj) ) face = 4
      if( (ks .eq. ke).and.(ks .eq.  1) ) face = 1
      if( (ks .eq. ke).and.(ks .eq. ck) ) face = 6
    end if

    ! Store boundary condition
    if(face .ne. 0) then

      found = .false.
      do r=1,n_bnd
        if( Grid % region % name(r) .eq.   &
            Dom % ranges(n) % name ) found = .true.
      end do
      if( .not. found) then
        n_bnd = n_bnd + 1
        Grid % region % name(n_bnd) = Dom % ranges(n) % name
      end if

      do i=is,ie
        do j=js,je
          do k=ks,ke
            c = Dom % blocks(b) % n_cells + (k-1)*ci*cj + (j-1)*ci + i
            Grid % cells_c(face,c) = -n_bnd
          end do
        end do
      end do

    end if

  end do  !  n_ranges

  ! Store the number of boundary conditions
  Grid % n_bnd_regions = n_bnd

  call Grid % Print_Regions_List()

  end subroutine
