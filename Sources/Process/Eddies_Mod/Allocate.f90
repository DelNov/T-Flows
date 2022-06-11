!==============================================================================!
  subroutine Eddies_Mod_Allocate(eddies,   &
                                 n_edd,    &
                                 max_r,    &
                                 intns,    &
                                 Flow, bnd_cond_name)
!------------------------------------------------------------------------------!
!   Allocate memory for eddies                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
  type(Field_Type),  target :: Flow
  integer                   :: n_edd   ! number of eddies
  real                      :: max_r   ! maximum eddy radius
  real                      :: intns   ! eddy intensity
  character(len=*)          :: bnd_cond_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: e
!==============================================================================!

  ! Store pointers to grid, flow and the name of boundary condition
  eddies % pnt_flow => Flow
  eddies % pnt_grid => Flow % pnt_grid
  eddies % bc_name  =  bnd_cond_name

  !---------------------!
  !   Allocate memory   !
  !---------------------!

  ! Store number of eddies ...
  eddies % n_eddies   = n_edd
  eddies % max_radius = max_r
  eddies % intensity  = intns

  ! ... and allocate memory for all of them
  allocate(eddies % eddy(eddies % n_eddies))

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!
  call Eddies_Mod_Gather_Bnd_Cells(eddies)

  !-------------------------------!
  !   Place all eddies randomly   !
  !-------------------------------!
  do e = 1, eddies % n_eddies
    call Eddies_Mod_Place_Eddy(eddies, e)
  end do

  end subroutine
