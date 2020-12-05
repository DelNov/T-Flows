!==============================================================================!
  program Convert
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: grid, dual     ! grid to be converted and its dual
  character(SL)   :: file_name
  character(SL)   :: file_format    ! 'UNKNOWN', 'FLUENT', 'GAMBIT', 'GMSH'
  integer         :: l, p
!==============================================================================!

  call Logo_Con

  print *, '#================================================================'
  print *, '# Enter the name of the mesh file you are importing (with ext.):'
  print *, '#----------------------------------------------------------------'
  read(*,*) file_name

  !-----------------------------------------------!
  !                                               !
  !   Make an educated guess of the file format   !
  !                                               !
  !-----------------------------------------------!
  call Guess_Format(file_name, file_format)

  !-------------------------------------------!
  !                                           !
  !   Extract the problem name and store it   !
  !                                           !
  !-------------------------------------------!
  l = len_trim(file_name)
  p = index(file_name(1:l), '.', back=.true.)

  problem_name(1) = file_name(1:p-1)

  grid % name = problem_name(1)
  call To_Upper_Case(grid % name)

  !----------------------------------------!
  !                                        !
  !   Read the file and start conversion   !
  !                                        !
  !----------------------------------------!
  if(file_format .eq. 'FLUENT') then
    call Load_Fluent(grid, file_name)
  end if
  if(file_format .eq. 'GAMBIT') then
    call Load_Gambit(grid, file_name)
  end if
  if(file_format .eq. 'GMSH') then
    call Load_Gmsh(grid, file_name)
    call Find_Parents(grid)
  end if

  ! For Gambit and Gmsh grids, no face information is stored
  if(file_format .eq. 'GAMBIT' .or. file_format .eq. 'GMSH') then
    call Grid_Topology(grid)
    call Find_Faces   (grid)
  end if
  call Calculate_Geometry(grid)

  ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
  ! ... faster if it is called after Grid_Mod_Sort_Faces_Smart
  call Grid_Mod_Sort_Cells_Smart       (grid)
  call Grid_Mod_Sort_Faces_Smart       (grid)
  call Grid_Mod_Calculate_Wall_Distance(grid)

  call Grid_Mod_Initialize_New_Numbers(grid)

  !-------------------------------!
  !                               !
  !   Save files for processing   !
  !                               !
  !-------------------------------!
  call Grid_Mod_Save_Cfn(grid, 0,             &
                         grid % n_nodes,      &
                         grid % n_cells,      &
                         grid % n_faces,      &
                         grid % n_shadows,    &
                         grid % n_bnd_cells)

  call Grid_Mod_Save_Geo(grid, 0)

  !-----------------------------------------------------!
  !                                                     !
  !   Save grid for visualisation and post-processing   !
  !                                                     !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)

  ! Create 1D file (used for channel or pipe flow)
  call Probe_1d_Nodes(grid)

  call Create_Dual(grid, dual)   ! It should be in Grid_Mod

  end program
