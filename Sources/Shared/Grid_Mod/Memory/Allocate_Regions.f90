!==============================================================================!
  subroutine Allocate_Regions(Grid, nbr)
!------------------------------------------------------------------------------!
!>  Allocates memory for boundary regions.
!------------------------------------------------------------------------------!
!   Regions.  There are some boundary regions for boundary conditions. These   !
!   are specified by a user in one file or another, depending if you are in    !
!   Generate, Convert or others.  They are stored in Grid % n_bnd_regions.     !
!                                                                              !
!   In addition to those, there is an extra region for internal cells or       !
!   faces.  It is stored in Grid % n_regions, and Grid % n_regions is equal    !
!   to Grid % n_bnd_regions + 1.                                               !
!                                                                              !
!   In addition, there are also buffer cells in Process, which would be        !
!   stored in Grid % n_regions + 1 or Grid % n_bnd_regions + 2.                !
!                                                                              !
!   On top of them all, there are three regions reserved for periodic x, y     !
!   and z directions, used for copy boundaries and they are stored in          !
!   Grid % n_regions + 2, Grid % n_regions + 3 and Grid % n_regions + 4.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: nbr   !! number of boundary regions
!==============================================================================!

  ! Store the number of boundary regions for the grid
  Grid % n_bnd_regions = nbr

  ! Number of regions, at least the "visible" ones ...
  ! ... is bigger by one, by the "internal" regaion
  Grid % n_regions = Grid % n_bnd_regions + 1

  Grid % per_x_reg = Grid % n_bnd_regions + 3   ! for periodic_x region
  Grid % per_y_reg = Grid % n_bnd_regions + 4   ! for periodic_y region
  Grid % per_z_reg = Grid % n_bnd_regions + 5   ! for periodic_z region

  ! Allocate memory for regions' names
  allocate(Grid % region % name(Grid % n_regions+5))
  Grid % region % name(1:Grid % n_bnd_regions)     = 'UNDEFINED'
  Grid % region % name  (Grid % n_bnd_regions + 1) = 'INSIDE'
  Grid % region % name  (Grid % n_bnd_regions + 2) = 'BUFFER'
  Grid % region % name  (Grid % n_bnd_regions + 3) = 'PERIODIC_X'
  Grid % region % name  (Grid % n_bnd_regions + 4) = 'PERIODIC_Y'
  Grid % region % name  (Grid % n_bnd_regions + 5) = 'PERIODIC_Z'

  ! Allocate memory for regions' types
  allocate(Grid % region % type(Grid % n_regions+5))
  Grid % region % type(1:Grid % n_bnd_regions)     = UNDEFINED
  Grid % region % type  (Grid % n_bnd_regions + 1) = INSIDE
  Grid % region % type  (Grid % n_bnd_regions + 2) = BUFFER
  Grid % region % type  (Grid % n_bnd_regions + 3) = PERIODIC_X
  Grid % region % type  (Grid % n_bnd_regions + 4) = PERIODIC_Y
  Grid % region % type  (Grid % n_bnd_regions + 5) = PERIODIC_Z

  ! Memory for cell and face ranges
  allocate(Grid % region % f_cell(Grid % n_regions+5))
  allocate(Grid % region % l_cell(Grid % n_regions+5))
  allocate(Grid % region % f_face(Grid % n_regions+5))
  allocate(Grid % region % l_face(Grid % n_regions+5))
  Grid % region % f_cell(:) = -1
  Grid % region % l_cell(:) = -HUGE_INT
  Grid % region % f_face(:) =  0
  Grid % region % l_face(:) = -1

  end subroutine
