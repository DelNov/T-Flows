!==============================================================================!
  program Convert
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod,      only: problem_name
  use Grid_Mod,      only: Grid_Type,                         &
                           Grid_Mod_Sort_Cells_Smart,         &
                           Grid_Mod_Sort_Faces_Smart,         &
                           Grid_Mod_Calculate_Wall_Distance,  &
                           Grid_Mod_Save_Cns,                 &
                           Grid_Mod_Save_Geo
  use Save_Grid_Mod, only: Save_Vtu_Cells,  &
                           Save_Vtu_Faces
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type)   :: grid     ! grid to be converted
  integer           :: c, n, s, l
  character(len=80) :: file_name, file_name_up, ext_up
!==============================================================================!

  call Logo_Con

  print *, '#==============================================================='
  print *, '# Enter the name of the mesh file you are importing (with ext.):'
  print *, '#---------------------------------------------------------------'
  read(*,*) file_name

  !-----------------------------------------------------!
  !   Analyze the extension to figure the file format   !
  !-----------------------------------------------------!
  file_name_up = file_name
  call To_Upper_Case(file_name_up)

  l = len_trim(file_name)
  print *, '#==================================' // &
           '==================================='
  if( file_name_up(l-2:l) .eq. 'NEU' ) then
    print *, '# Based on the extension, you are' // &
             ' reading Gambit''s neutral file format'
    problem_name(1) = file_name(1:l-4)
    ext_up = file_name_up(l-2:l)
  else if( file_name_up(l-2:l) .eq. 'MSH' ) then
    print *, '# Based on the extension, you are' // &
             ' reading GMSH file format'
    problem_name(1) = file_name(1:l-4)
    ext_up = file_name_up(l-2:l)
  else
    print *, '# Unrecognized input file format; exiting!'
    print *, '#----------------------------------' // &
             '-----------------------------------'
    stop
  end if
  print *, '#----------------------------------' // &
           '-----------------------------------'

  grid % name = problem_name(1)
  call To_Upper_Case(grid % name)

  !----------------------------------------!
  !   Read the file and start conversion   !
  !----------------------------------------!
  if (ext_up .eq. 'NEU') then
    call Load_Neu(grid)
  end if
  if (ext_up .eq. 'MSH') then
    call Load_Msh(grid)
    call Find_Parents(grid)
  end if

  call Grid_Topology     (grid)
  call Find_Faces        (grid)
  call Calculate_Geometry(grid)

  ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
  ! ... faster if it is called after Grid_Mod_Sort_Faces_Smart

  call Grid_Mod_Sort_Cells_Smart       (grid)
  call Grid_Mod_Sort_Faces_Smart       (grid)
  call Grid_Mod_Calculate_Wall_Distance(grid)

  ! Prepare for saving
  do n = 1, grid % n_nodes
    grid % new_n(n) = n
  end do
  do c = -grid % n_bnd_cells, grid % n_cells
    grid % new_c(c) = c
    grid % old_c(c) = c
  end do
  do s = 1, grid % n_faces + grid % n_shadows
    grid % new_f(s) = s
    grid % old_f(s) = s
  end do

  ! Decompose/coarsen the grid with METIS
  ! Not quite ready yet: call Grid_Mod_Coarsen(grid)

  !-------------------------------!
  !   Save files for processing   !
  !-------------------------------!
  call Grid_Mod_Save_Cns(grid, 0,             &
                         grid % n_nodes,      &
                         grid % n_cells,      &
                         grid % n_faces,      &
                         grid % n_shadows,    &
                         grid % n_bnd_cells)

  call Grid_Mod_Save_Geo(grid, 0)

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

 ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)

  ! Create 1D file (used for channel or pipe flow)
  call Probe_1d_Nodes(grid)

  end program
