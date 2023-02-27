!==============================================================================!
  subroutine Determine_Regions_Ranges(Grid)
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
  integer :: c, c1, c2, reg, s, siz
!==============================================================================!

  !-------------------!
  !   Cells' ranges   !
  !-------------------!

  ! Set non-realizable ranges
  Grid % region % f_cell(:) = -1
  Grid % region % l_cell(:) = -HUGE_INT

  ! Browse forward and backward to find first and last cell for each range
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      reg = Grid % region % at_cell(c)
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
    end if
  end do

  do c = -1, -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      reg = Grid % region % at_cell(c)
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  ! Inside cells
  reg = Grid % n_regions
  Grid % region % f_cell(reg) = Grid % n_cells
  Grid % region % l_cell(reg) = 1
  do c = 1, Grid % n_cells
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  ! Buffer cells
  reg = Grid % n_regions + 1
  Grid % region % f_cell(reg) = Grid % n_cells + 1
  Grid % region % l_cell(reg) = 1
  do c = 1, Grid % n_cells
    if(Grid % Comm % cell_proc(c) .ne. this_proc) then
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  if(DEBUG) then
    write(1000+this_proc, '(a,a)')  ' # Cell ranges from ', Program_Name
    do reg = All_Regions()
      siz = Grid % region % l_cell(reg) - Grid % region % f_cell(reg) + 1
      write(1000+this_proc,'(a,i3,i15,i15,i15,a,a)') ' # Region: ', reg, &
                                           Grid % region % f_cell(reg),  &
                                           Grid % region % l_cell(reg),  &
                                           max(siz, 0), '  ',            &
                                           trim(Grid % region % name(reg))
    end do
  end if

  !-------------------!
  !   Faces' ranges   !
  !-------------------!

  ! Set non-realizable ranges
  Grid % region % f_face(:) =  0
  Grid % region % l_face(:) = -1

  ! Browse through colors and faces to set faces' ranges
  do reg = Boundary_Regions()
    if(Grid % region % f_cell(reg) < Grid % region % l_cell(reg)) then
      do s = 1, Grid % n_faces
        c = Grid % faces_c(2, s)
        if(Grid % region % f_cell(reg) .eq. c) then
          Grid % region % f_face(reg) = s
        end if
        if(Grid % region % l_cell(reg) .eq. c) then
          Grid % region % l_face(reg) = s
        end if
      end do
    end if
  end do

  ! Inside faces
  reg = Grid % n_regions
  Grid % region % f_face(reg) = Grid % n_faces
  Grid % region % l_face(reg) = 1
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 > 0) then  ! limit to inside faces
      if(Grid % Comm % cell_proc(c1) .eq. this_proc .or.  &
         Grid % Comm % cell_proc(c2) .eq. this_proc) then
        if(s <= Grid % region % f_face(reg)) then
          Grid % region % f_face(reg) = s
        end if
        if(s >= Grid % region % l_face(reg)) then
          Grid % region % l_face(reg) = s
        end if
        Assert(Grid % Comm % cell_proc(c1) .eq. this_proc)  ! this should hold
      end if
    end if  ! c2 > 0
  end do

  do s = Faces_In_Domain()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    Assert(c2 > 0)
  end do

  if(DEBUG) then
    write(2000+this_proc, '(a,a)')  ' # Face ranges from ', Program_Name
    do reg = All_Regions()
      siz = Grid % region % l_face(reg) - Grid % region % f_face(reg) + 1
      write(2000+this_proc,'(a,i3,i15,i15,i15,a,a)') ' # Region: ', reg, &
                                           Grid % region % f_face(reg),  &
                                           Grid % region % l_face(reg),  &
                                           max(siz, 0), '  ',            &
                                           trim(Grid % region % name(reg))
    end do
  end if

  !-----------!
  !   Check   !
  !-----------!
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c = Grid % faces_c(2, s)
      Assert(Grid % region % at_cell(c) == reg)
    end do
  end do

  end subroutine
