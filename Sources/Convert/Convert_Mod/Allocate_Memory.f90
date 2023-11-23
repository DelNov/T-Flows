!==============================================================================!
  subroutine Allocate_Memory(Convert, Grid)
!------------------------------------------------------------------------------!
!>  The subroutine is responsible for dynamically allocating memory for
!>  different elements of the grid such as nodes, cells, and faces, for the
!>  Convert sub-program.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: F = 3   ! workaround to allocate more memory for bnds
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Allocate memory
  ! =--> carefull: there is no checking!
  call Grid % Allocate_Nodes(Grid % n_nodes)
  call Grid % Allocate_Cells(Grid % n_cells,   Grid % n_bnd_cells * F)

  ! For Gmsh's and Gambit's file formats, number of faces is zero
  if(Grid % n_faces .eq. 0) then
    call Grid % Allocate_Faces(Grid % n_cells*5, 0)

  ! For Fluent's file format, number of faces are known
  else
    call Grid % Allocate_Faces(Grid % n_faces, 0)
  end if

  allocate(Grid % region % at_cell(-Grid % n_bnd_cells * F:-1))
  Grid % region % at_cell = 0

  ! (Dirty) trick to allocate additional memory for cities
  Grid % n_bnd_cells = Grid % n_bnd_cells / F

  end subroutine
