!==============================================================================!
  subroutine Determine_Regions_Ranges(Grid)
!------------------------------------------------------------------------------!
!   Allocates memory and finds the range (first and last boundary cell)        !
!   for each of the boundary condition region.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, reg, s, siz
  integer :: last_face_only_inside, first_face_in_buffers
!==============================================================================!

  !-------------------!
  !                   !
  !   Cells' ranges   !
  !                   !
  !-------------------!

  !-------------------------------!
  !   Set non-realizable ranges   !
  !-------------------------------!
  Grid % region % f_cell(:) = -1
  Grid % region % l_cell(:) = -HUGE_INT

  !----------------------------------------------------------------------------!
  !   Browse forward and backward to find first and last cell for each range   !
  !----------------------------------------------------------------------------!

  ! Forward
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      reg = Grid % region % at_cell(c)
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
    end if
  end do

  ! Backward
  do c = -1, -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      reg = Grid % region % at_cell(c)
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  !------------------!
  !   Inside cells   !
  !------------------!
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

  !------------------!
  !   Buffer cells   !
  !- - - - - - - - - +--------------------------------------------------!
  !   Note: leave Grid % region % l_cell(reg) at Grid % n_cells since   !
  !   it makes macro Cells_In_Domain_And_Buffers() work in sequential   !
  !---------------------------------------------------------------------!
  reg = Grid % n_regions + 1
  Grid % region % f_cell(reg) = Grid % n_cells + 1
  Grid % region % l_cell(reg) = Grid % n_cells
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
    write(1000+this_proc, '(a,a)')  ' # Cell ranges from ', PROGRAM_NAME
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
  !                   !
  !   Faces' ranges   !
  !                   !
  !-------------------!

  !-------------------------------!
  !   Set non-realizable ranges   !
  !-------------------------------!
  Grid % region % f_face(:) =  0
  Grid % region % l_face(:) = -1

  !----------------------------------------------------------!
  !   Browse through colors and faces to set faces' ranges   !
  !----------------------------------------------------------!
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

  !------------------!
  !   Inside faces   !
  !------------------!
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

  !-------------!
  !   Check 1   !
  !-------------!
  do s = Faces_In_Domain()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    Assert(c2 > 0)
  end do

  !-------------!
  !   Check 2   !
  !-------------!
  last_face_only_inside = 0
  first_face_in_buffers = HUGE_INT  ! first face in buffers
  do s = Faces_In_Domain()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(Grid % Comm % cell_proc(c1) .eq. this_proc .and.  &
       Grid % Comm % cell_proc(c2) .eq. this_proc) then
      last_face_only_inside = max(last_face_only_inside, s)
    end if
    if(Grid % Comm % cell_proc(c1) .eq. this_proc .and.  &
       Grid % Comm % cell_proc(c2) .ne. this_proc) then
      first_face_in_buffers = min(first_face_in_buffers, s)
    end if
  end do
  Assert(last_face_only_inside < first_face_in_buffers)

  if(DEBUG) then
    write(2000+this_proc, '(a,a)')  ' # Face ranges from ', PROGRAM_NAME
    do reg = All_Regions()
      siz = Grid % region % l_face(reg) - Grid % region % f_face(reg) + 1
      write(2000+this_proc,'(a,i3,i15,i15,i15,a,a)') ' # Region: ', reg, &
                                           Grid % region % f_face(reg),  &
                                           Grid % region % l_face(reg),  &
                                           max(siz, 0), '  ',            &
                                           trim(Grid % region % name(reg))
    end do
  end if

  !-------------!
  !   Check 3   !
  !-------------!
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c = Grid % faces_c(2, s)
      Assert(Grid % region % at_cell(c) == reg)
    end do
  end do

  end subroutine
