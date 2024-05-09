!==============================================================================!
  subroutine Save_Fluent(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to save grid files in the Fluent's file format
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File opening and initial setup: The subroutine opens the Fluent file in  !
!     binary mode, assuming that the grid does not contain polyhedral cells    !
!     initially.                                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being saved
!----------------------------------[Calling]-----------------------------------!
# ifdef __NVCOMPILER_LLVM__
  integer :: ftell  ! needed for Nvidia compiler
# endif
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MIXED_ZONE = 0
  integer, parameter :: CELL_TRI   = 1
  integer, parameter :: CELL_TETRA = 2
  integer, parameter :: CELL_QUAD  = 3
  integer, parameter :: CELL_HEXA  = 4
  integer, parameter :: CELL_PYRA  = 5
  integer, parameter :: CELL_WEDGE = 6
  integer, parameter :: CELL_POLY  = 7
  integer, parameter :: FACE_TRI   = 3
  integer, parameter :: FACE_QUAD  = 4
  integer, parameter :: FACE_POLY  = 5  ! an educated guess, but seems good
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: file_name
  integer       :: fu, n, c, i_nod, s, c1, c2, reg, reg_type, reg_zone
!==============================================================================!

  call Profiler % Start('Save_Fluent')

  !------------------------------------------!
  !                                          !
  !   Set the name of the file and open it   !
  !                                          !
  !------------------------------------------!
  file_name = trim(Grid % name)
  call String % To_Lower_Case(file_name)
  call File % Set_Name(file_name, extension=".cas")
  call File % Open_For_Writing_Ascii(file_name, fu)

  !------------------!
  !                  !
  !   Write header   !
  !                  !
  !------------------!

  !------------------------------------------------!
  !   Dimensions; it is a three-dimensional grid   !
  !------------------------------------------------!
  write(fu, '(a)') "(2 3)"

  !--------------------------------------------!
  !   Info on the software creating the file   !
  !--------------------------------------------!
  write(fu, '(a)') "(0 ""Exported from T-Flows"")"
  write(fu, '(a)') ""

  !--------------------------------------------!
  !   Write number of nodes, cells and faces   !
  !--------------------------------------------!
  write(fu, '(a, z9, a)')  "(10 (0 1 ", Grid % n_nodes, " 0))"
  write(fu, '(a, z9, a)')  "(12 (0 1 ", Grid % n_cells, " 0))"
  write(fu, '(a, z9, a)')  "(13 (0 1 ", Grid % n_faces, " 0))"

  !-----------------!
  !                 !
  !   Write nodes   !
  !                 !
  !-----------------!
  write(fu, '(a, z9, a)')  "(10 (1 1 ", Grid % n_nodes, " 1) ("
  do n = 1, Grid % n_nodes
    write(fu, '(3es15.6)')  Grid % xn(n), Grid % yn(n), Grid % zn(n)
  end do
  write(fu, '(a)')  "))"

  !-----------------!
  !                 !
  !   Write cells   !
  !                 !
  !-----------------!
  write(fu, '(a, 3z9, a)') "(12 (1 1 ", Grid % n_cells, 1, MIXED_ZONE, ") ("

  do c = 1, Grid % n_cells
    select case (Grid % cells_n_nodes(c))
      case (4)
        write(fu, '(z2)', advance="no")  CELL_TETRA
      case (5)
        write(fu, '(z2)', advance="no")  CELL_PYRA
      case (6)
        write(fu, '(z2)', advance="no")  CELL_WEDGE
      case (8)
        write(fu, '(z2)', advance="no")  CELL_HEXA
      case default
        PRINT *, 'ERROR!'
        STOP
    end select
  end do
  write(fu, '(a)')  "))"

  !-----------------!
  !                 !
  !   Write faces   !
  !                 !
  !-----------------!
  do reg = 1, Grid % n_bnd_regions + 1

    ! Specify region type
    reg_type = 3                              ! default is wall
    if(reg .eq. Grid % n_bnd_regions + 1)  &  ! correct for interior
      reg_type = 2

    ! Specify region zone
    reg_zone = reg + 2                        ! default are boundary regions
    if(reg .eq. Grid % n_bnd_regions + 1)  &  ! correct for interiors
      reg_zone = 2

    ! Leading line of the new region
    write(fu, '(a, z0, z9, z9, z2, z2, a)')     &
               "(13 (",                        &
               reg_zone,                       &
               Grid % region % f_face(reg),    &
               Grid % region % l_face(reg),    &
               reg_type,                       &
               MIXED_ZONE,                     &
               ") ("

    !------------------------------------------!
    !   Browse through faces of that regions   !
    !------------------------------------------!
    do s = Grid % region % f_face(reg), Grid % region % l_face(reg)

      ! Write the number of nodes for this face
      write(fu, '(z9, 1x)', advance="no")  Grid % faces_n_nodes(c)

      ! Write all nodes for the face
      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        write(fu, '(z9, 1x)', advance="no")  n
      end do

      ! Write cells surroundig the face
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      c2 = max(0, c2)            ! if it is a boundary cell, set to zero

      write(fu, '(2z9)', advance="no")  c1, c2

      ! End the line for face s
      write(fu, '()')
    end do
    write(fu, '(a)')  "))"
  end do

  !--------------------!
  !                    !
  !   Close the file   !
  !                    !
  !--------------------!
  close(fu)

  call Profiler % Stop('Save_Fluent')

  end subroutine
