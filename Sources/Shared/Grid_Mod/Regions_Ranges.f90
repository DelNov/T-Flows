!==============================================================================!
  subroutine Regions_Ranges(Grid)
!------------------------------------------------------------------------------!
!   Allocates memory and finds the range (first and last boundary cell)        !
!   for each of the boundary condition region.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: c2, reg, s
!==============================================================================!

  ! Allocate memory
  allocate(Grid % region % f_cell(0:Grid % n_regions))
  allocate(Grid % region % l_cell(0:Grid % n_regions))

  !-------------------!
  !   Cells' ranges   !
  !-------------------!

  ! Set non-realizable ranges
  Grid % region % f_cell(:) = -1
  Grid % region % l_cell(:) = -HUGE_INT

  ! Browse forward and backward to find first and last cell for each range
  do c2 = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c2) .eq. this_proc) then
      reg = Grid % region % at_cell(c2)
      if(c2 < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c2
      end if
    end if
  end do

  do c2 = -1, -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c2) .eq. this_proc) then
      reg = Grid % region % at_cell(c2)
      if(c2 > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c2
      end if
    end if
  end do

  if(DEBUG) then
    write(1000, '(a)')  ' # Cell ranges'
    do reg = 1, Grid % n_regions
      write(1000+this_proc,'(a,i3,i9,i9)') ' # Region: ', reg,   &
                                   Grid % region % f_cell(reg),  &
                                   Grid % region % l_cell(reg)
    end do
  end if

  !-------------------!
  !   Faces' ranges   !
  !-------------------!

  ! Allocate memory
  allocate(Grid % region % f_face(0:Grid % n_regions))
  allocate(Grid % region % l_face(0:Grid % n_regions))

  ! Set non-realizable ranges
  Grid % region % f_face(:) = -1
  Grid % region % l_face(:) = -HUGE_INT

  ! Browse through colors and faces to set faces' ranges
  do reg = 1, Grid % n_regions
    do s = 1, Grid % n_faces
      c2 = Grid % faces_c(2, s)
      if(Grid % region % f_cell(reg) .eq. c2) then
        Grid % region % f_face(reg) = s
      end if
      if(Grid % region % l_cell(reg) .eq. c2) then
        Grid % region % l_face(reg) = s
      end if
    end do
  end do

  if(DEBUG) then
    write(2000, '(a)')  ' # Face ranges'
    do reg = 1, Grid % n_regions
      write(2000+this_proc,'(a,i3,i9,i9)') ' # Region: ', reg,   &
                                   Grid % region % f_face(reg),  &
                                   Grid % region % l_face(reg)
    end do
  end if

  !-----------!
  !   Check   !
  !-----------!
  do reg = 1, Grid % n_regions
    Browse_Faces_In_Color(s, reg)
      c2 = Grid % faces_c(2, s)
      Assert(Grid % region % at_cell(c2) == reg)
    End_Browse
  end do

  end subroutine
