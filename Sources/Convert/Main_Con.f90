!==============================================================================!
  program Convert
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Grid_Mod, only: Grid_Type,                        &
                      Grid_Mod_Sort_Faces_Smart,        &
                      Grid_Mod_Calculate_Wall_Distance
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type)   :: grid     ! grid to be converted
  integer           :: c, n, s, l, color
  character(len=80) :: file_name, file_name_up, extension
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
    problem_name = file_name(1:l-4)
    extension = file_name_up(l-2:l)
  else if( file_name_up(l-3:l) .eq. 'CGNS' ) then
    print *, '# Based on the extension, you are' // &
             ' reading CGNS file format'
    problem_name = file_name(1:l-5)
    extension = file_name_up(l-3:l)
  else
    print *, '# Unrecognized input file format; exiting!'
    stop
  end if
  print *, '#----------------------------------' // &
           '-----------------------------------'

  !----------------------------------------!
  !   Read the file and start conversion   !
  !----------------------------------------!
  if (extension .eq. 'NEU') then
    call Load_Neu (grid)
  end if
  if (extension .eq. 'CGNS') then
    call Load_Cgns(grid)
  end if

  call Grid_Topology     (grid)
  call Find_Faces        (grid)
  call Calculate_Geometry(grid)
  call Connect_Domains   (grid)

  ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
  ! ... faster if it is called after Grid_Mod_Sort_Faces_Smart
  call Grid_Mod_Sort_Faces_Smart       (grid)
  call Grid_Mod_Calculate_Wall_Distance(grid)

  ! Prepare for saving
  do n = 1, grid % n_nodes
    grid % new_n(n) = n
  end do
  do c = -grid % n_bnd_cells, grid % n_cells
    grid % new_c(c) = c
  end do
  do s = 1, grid % n_faces
    grid % new_f(s) = s
  end do

  ! Decompose/coarsen the grid with METIS
  ! Not yet implemented: call Grid_Mod_Coarsen(grid)

  !-------------------------------!
  !   Save files for processing   !
  !-------------------------------!
  call Save_Cns_Geo(grid, 0,             &
                    grid % n_nodes,      &
                    grid % n_cells,      &
                    grid % n_faces,      &
                    grid % n_bnd_cells,  &
                    0, 0)

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

 ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid, 0)

  ! Save links for checking
  call Save_Vtu_Links(grid, 0,             &
                      grid % n_nodes,      &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0)

  ! Create 1D file (used for channel or pipe flow)
  call Probe_1d_Nodes(grid)

  end program
