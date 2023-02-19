!==============================================================================!
  subroutine Bnd_Cond_Ranges(Grid)
!------------------------------------------------------------------------------!
!   Allocates memory and finds the range (first and last boundary cell)        !
!   for each of the boundary condition colors.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: c2, color, s
!==============================================================================!

  ! Allocate memory
  allocate(Grid % bnd_cond % color_s_cell(0:Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % color_e_cell(0:Grid % n_bnd_cond))

  !-------------------!
  !   Cells' ranges   !
  !-------------------!

  ! Set non-realizable ranges
  Grid % bnd_cond % color_s_cell(:) = -1
  Grid % bnd_cond % color_e_cell(:) = -HUGE_INT

  ! Browse forward and backward to find first and last cell for each range
  do c2 = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c2) .eq. this_proc) then
      color = Grid % bnd_cond % color(c2)
      if(c2 < Grid % bnd_cond % color_s_cell(color)) then
        Grid % bnd_cond % color_s_cell(color) = c2
      end if
    end if
  end do

  do c2 = -1, -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c2) .eq. this_proc) then
      color = Grid % bnd_cond % color(c2)
      if(c2 > Grid % bnd_cond % color_e_cell(color)) then
        Grid % bnd_cond % color_e_cell(color) = c2
      end if
    end if
  end do

  if(DEBUG) then
    write(1000, '(a)')  ' # Cell ranges'
    do color = 1, Grid % n_bnd_cond
      write(1000+this_proc,'(a,i3,i9,i9)') ' # Region: ', color,  &
                           Grid % bnd_cond % color_s_cell(color),  &
                           Grid % bnd_cond % color_e_cell(color)
    end do
  end if

  !-------------------!
  !   Faces' ranges   !
  !-------------------!

  ! Allocate memory
  allocate(Grid % bnd_cond % color_s_face(0:Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % color_e_face(0:Grid % n_bnd_cond))

  ! Set non-realizable ranges
  Grid % bnd_cond % color_s_face(:) = -1
  Grid % bnd_cond % color_e_face(:) = -HUGE_INT

  ! Browse through colors and faces to set faces' ranges
  do color = 1, Grid % n_bnd_cond
    do s = 1, Grid % n_faces
      c2 = Grid % faces_c(2, s)
      if(Grid % bnd_cond % color_s_cell(color) .eq. c2) then
        Grid % bnd_cond % color_s_face(color) = s
      end if
      if(Grid % bnd_cond % color_e_cell(color) .eq. c2) then
        Grid % bnd_cond % color_e_face(color) = s
      end if
    end do
  end do

  if(DEBUG) then
    write(2000, '(a)')  ' # Face ranges'
    do color = 1, Grid % n_bnd_cond
      write(2000+this_proc,'(a,i3,i9,i9)') ' # Region: ', color,  &
                           Grid % bnd_cond % color_s_face(color),  &
                           Grid % bnd_cond % color_e_face(color)
    end do
  end if

  !-----------!
  !   Check   !
  !-----------!
  do color = 1, Grid % n_bnd_cond
    Browse_Faces_In_Color(s, color)
      c2 = Grid % faces_c(2, s)
      Assert(Grid % bnd_cond % color(c2) == color)
    End_Browse
  end do

  end subroutine
