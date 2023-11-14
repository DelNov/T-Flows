!==============================================================================!
  subroutine Determine_Regions_Ranges(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine identifies the ranges for each of the boundary condition
!>  regions within the computational grid. It is a key function in setting up
!>  the grid for simulations, particularly useful in Process, where it allows
!>  for efficient navigation through different regions of the grid.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization of cell ranges:                                           !
!     - Sets initial non-realizable ranges for cells within each boundary      !
!       region.                                                                !
!   * Forward and Backward Search:                                             !
!     - Iterates through the grid cells in both forward and backward direction !
!       to find the first and last cell for each boundary condition region.    !
!   * Inside Cells:                                                            !
!     - Determines the range of cells inside the domain (excluding buffer      !
!       cells).                                                                !
!   * Buffer Cells:                                                            !
!     - Identifies the range of buffer cells, which are important for parallel !
!       processing and communication between subdomains.                       !
!   * Debugging:                                                               !
!     - If the DEBUG mode is enabled, outputs the range of cells for each      !
!       region to assist in the debugging process.                             !
!   * Face Ranges:                                                             !
!     - Similar to cells, it sets ranges for faces within each region. This    !
!       includes boundary faces, inside faces, and faces in the buffer.        !
!   * Verification:                                                            !
!     - Performs checks to ensure the integrity and correctness of the ranges  !
!       set for both cells and faces in the grid.                              !
!   * Utilization of Macros:                                                   !
!     - Once this subroutine is exexcuted, the macros defined in in Browse.h90 !
!       can be utilizied for safer browsing and iteration through regions.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! computational grid
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
    if(Cell_In_This_Proc(c)) then
      reg = Grid % region % at_cell(c)
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
    end if
  end do

  ! Backward
  do c = -1, -Grid % n_bnd_cells, -1
    if(Cell_In_This_Proc(c)) then
      reg = Grid % region % at_cell(c)
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  !------------------!
  !   Inside cells   !
  !- - - - - - - - - +--------------------------------------!
  !   Note that this does not entail cells in the buffers   !
  !---------------------------------------------------------!
  reg = Grid % n_regions
  Grid % region % f_cell(reg) = Grid % n_cells
  Grid % region % l_cell(reg) = 1
  do c = 1, Grid % n_cells
    if(Cell_In_This_Proc(c)) then
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
    if(.not. Cell_In_This_Proc(c)) then
      if(c < Grid % region % f_cell(reg)) then
        Grid % region % f_cell(reg) = c
      end if
      if(c > Grid % region % l_cell(reg)) then
        Grid % region % l_cell(reg) = c
      end if
    end if
  end do

  if(DEBUG) then
    write(1000+This_Proc(), '(a,a)')  ' # Cell ranges from ', PROGRAM_NAME
    do reg = All_Regions()
      siz = Grid % region % l_cell(reg) - Grid % region % f_cell(reg) + 1
      write(1000+This_Proc(),'(a,i3,i15,i15,i15,a,a)') ' # Region: ', reg, &
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
  !- - - - - - - - - +--------------------------------------------------!
  !   Note: Unlike inside cells, this also holds faces in the buffers   !
  !---------------------------------------------------------------------!
  reg = Grid % n_regions
  Grid % region % f_face(reg) = Grid % n_faces
  Grid % region % l_face(reg) = 1
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 > 0) then  ! limit to inside faces
      if(Cell_In_This_Proc(c1) .or. Cell_In_This_Proc(c2)) then
        if(s <= Grid % region % f_face(reg)) then
          Grid % region % f_face(reg) = s
        end if
        if(s >= Grid % region % l_face(reg)) then
          Grid % region % l_face(reg) = s
        end if
        Assert(Cell_In_This_Proc(c1))  ! this should hold
      end if
    end if  ! c2 > 0
  end do

  !-------------!
  !   Check 1   !
  !-------------!
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    Assert(c2 > 0)
  end do

  !-------------!
  !   Check 2   !
  !-------------!
  last_face_only_inside = 0
  first_face_in_buffers = HUGE_INT  ! first face in buffers
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(Cell_In_This_Proc(c1) .and. Cell_In_This_Proc(c2)) then
      last_face_only_inside = max(last_face_only_inside, s)
    end if
    if(Cell_In_This_Proc(c1) .and. .not. Cell_In_This_Proc(c2)) then
      first_face_in_buffers = min(first_face_in_buffers, s)
    end if
  end do
  Assert(last_face_only_inside < first_face_in_buffers)

  if(DEBUG) then
    write(2000+This_Proc(), '(a,a)')  ' # Face ranges from ', PROGRAM_NAME
    do reg = All_Regions()
      siz = Grid % region % l_face(reg) - Grid % region % f_face(reg) + 1
      write(2000+This_Proc(),'(a,i3,i15,i15,i15,a,a)') ' # Region: ', reg,    &
                                             Grid % region % f_face(reg),     &
                                             Grid % region % l_face(reg),     &
                                             max(siz, 0), '  ',               &
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
