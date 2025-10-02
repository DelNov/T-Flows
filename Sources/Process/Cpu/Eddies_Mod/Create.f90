!==============================================================================!
  subroutine Create_Eddies(Eddies,  &
                           n_edd,   &
                           max_r,   &
                           intns,   &
                           Flow,    &
                           bnd_cond_name)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to allocate memory and initialize eddies in a
!>  CFD simulation for creating synthetic turbulence at the inlet.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!   * Initialization and memory allocation:                                    !
!     - The subroutine first links the eddies structure to the flow and its    !
!       associated grid.                                                       !
!     - It then stores the name of the boundary condition and the number of    !
!       eddies, maximum eddy radius, and eddy intensity.                       !
!     - Memory is allocated for the specified number of eddies in the          !
!       eddies % eddy array.                                                   !
!   * Eddy placement:                                                          !
!     - The subroutine gathers information about boundary cells corresponding  !
!       to the specified boundary condition using the Gather_Bnd_Cells call.   !
!       This information is used to position eddies relative to these boundary !
!       cells.                                                                 !
!     - Each eddy is then placed randomly within the domain by calling         !
!       Eddies_Mod_Place_Eddy for each eddy. This involves setting the eddy's  !
!       initial position, radius, velocity, and other properties.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type),   target :: Eddies         !! parent class;
                                                 !! collection of eddies
  type(Field_Type),     target :: Flow           !! flow field on which
                                                 !! to impose eddies
  integer,          intent(in) :: n_edd          !! number of eddies
  real,             intent(in) :: max_r          !! maximum eddy radius
  real,             intent(in) :: intns          !! eddy intensity
  character(len=*), intent(in) :: bnd_cond_name  !! name of the boundary
                                                 !! condition to attach eddies
!-----------------------------------[Locals]-----------------------------------!
  integer :: e
!==============================================================================!

  ! Store pointers to grid, flow and the name of boundary condition
  Eddies % pnt_flow => Flow
  Eddies % pnt_grid => Flow % pnt_grid
  Eddies % bc_name  =  bnd_cond_name

  !---------------------!
  !   Allocate memory   !
  !---------------------!

  ! Store number of eddies ...
  Eddies % n_eddies   = n_edd
  Eddies % max_radius = max_r
  Eddies % intensity  = intns

  ! ... and allocate memory for all of them
  allocate(Eddies % eddy(Eddies % n_eddies))

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!
  call Eddies % Gather_Bnd_Cells()

  !-------------------------------!
  !   Place all eddies randomly   !
  !-------------------------------!
  do e = 1, Eddies % n_eddies
    call Eddies % Place_Eddy(e)
  end do

  end subroutine
