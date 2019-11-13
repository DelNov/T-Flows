!==============================================================================!
  subroutine Eddies_Mod_Allocate(eddies, n, max_r, flow, bnd_cond_name)
!------------------------------------------------------------------------------!
!   Allocate memory for eddies                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
  type(Field_Type),  target :: flow
  integer                   :: n       ! number of eddies
  real                      :: max_r   ! maximum eddy radius
  character(len=*)          :: bnd_cond_name
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: e, c, run
  real                     :: tmp
!==============================================================================!

  ! Store useful pointers
  eddies % pnt_flow => flow
  eddies % pnt_grid => flow % pnt_grid

  grid => flow % pnt_grid

  ! Store number of eddies ...
  eddies % n_eddies   = n
  eddies % bc_name    = bnd_cond_name
  eddies % max_radius = max_r

  ! ... and allocate memory for all of them
  allocate(eddies % eddy(eddies % n_eddies))

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!
  do run = 1, 2
    if(run .eq. 2) allocate(eddies % bnd_cell(eddies % n_bnd_cells))
    eddies % n_bnd_cells = 0
    do c = -grid % n_bnd_cells, -1
      if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then
        eddies % n_bnd_cells = eddies % n_bnd_cells + 1
        if(run .eq. 2) eddies % bnd_cell(eddies % n_bnd_cells) = c
      end if
    end do
  end do
  print '(a,a,i6)', ' # Number of boundary cells at ',  &
                    trim(eddies % bc_name),             &
                    eddies % n_bnd_cells

  !-------------------------------!
  !   Place all eddies randomly   !
  !-------------------------------!
  do e = 1, eddies % n_eddies
    call Eddies_Mod_Place_Eddy(eddies, e)
  end do

  end subroutine
