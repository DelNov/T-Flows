!==============================================================================!
  program Convert
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: grid(2)        ! grid to be converted and its dual
  character(SL)   :: answer
  character(SL)   :: file_name
  character(SL)   :: file_format    ! 'UNKNOWN', 'FLUENT', 'GAMBIT', 'GMSH'
  integer         :: l, p, g, n_grids
!==============================================================================!

  call Logo_Con

  print *, '#================================================================'
  print *, '# Enter the name of the grid file you are importing (with ext.):'
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

  grid(1) % name = problem_name(1)
  call To_Upper_Case(grid(1) % name)

  !----------------------------------------!
  !                                        !
  !   Read the file and start conversion   !
  !                                        !
  !----------------------------------------!
  if(file_format .eq. 'FLUENT') then
    call Load_Fluent(grid(1), file_name)
  end if
  if(file_format .eq. 'GAMBIT') then
    call Load_Gambit(grid(1), file_name)
  end if
  if(file_format .eq. 'GMSH') then
    call Load_Gmsh(grid(1), file_name)
    call Find_Parents(grid(1))
  end if

  ! For Gambit and Gmsh grids, no face information is stored
  if(file_format .eq. 'GAMBIT' .or. file_format .eq. 'GMSH') then
    call Grid_Topology(grid(1))
    call Find_Faces   (grid(1))
  end if

  !---------------------------------------------------!
  !                                                   !
  !   Decide if you are going for dual grid as well   !
  !                                                   !
  !---------------------------------------------------!
  print *, '#================================================='
  print *, '# Would you like to create a dual grid? (yes/no)'
  print *, '#-------------------------------------------------'
  call File_Mod_Read_Line(5)
  answer = line % tokens(1)
  call To_Upper_Case(answer)

  n_grids = 1
  if(answer .eq. 'YES') n_grids = 2

  !-------------------------------!
  !                               !
  !   Browse through both grids   !
  !                               !
  !-------------------------------!
  do g = 1, n_grids

    if(n_grids .eq. 2) then
      print *,              '#======================================'
      print *,              '#                                      '
      if(g .eq. 1) print *, '# Processing the first (primal) grid'
      if(g .eq. 2) print *, '# Processing the second (dual) grid'
      print *,              '#                                    '
      print *,              '#--------------------------------------'
      if(g .eq. 2) call Create_Dual(grid(1), grid(2))
    end if

    !--------------------------------------!
    !   Calculate geometrical quantities   !
    !--------------------------------------!
    call Calculate_Geometry(grid(g))

    ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
    ! ... faster if it is called after Grid_Mod_Sort_Faces_Smart
    call Grid_Mod_Sort_Cells_Smart       (grid(g))
    call Grid_Mod_Sort_Faces_Smart       (grid(g))
    call Grid_Mod_Calculate_Wall_Distance(grid(g))
    call Grid_Mod_Find_Cells_Faces       (grid(g))

    call Grid_Mod_Initialize_New_Numbers(grid(g))

    !-------------------------------!
    !   Save files for processing   !
    !-------------------------------!
    call Grid_Mod_Save_Cfn(grid(g), 0,             &
                           grid(g) % n_nodes,      &
                           grid(g) % n_cells,      &
                           grid(g) % n_faces,      &
                           grid(g) % n_shadows,    &
                           grid(g) % n_bnd_cells)

    call Grid_Mod_Save_Dim(grid(g), 0)

    !-----------------------------------------------------!
    !   Save grid for visualisation and post-processing   !
    !-----------------------------------------------------!

    ! Create output in vtu format
    call Save_Vtu_Cells(grid(g), 0,         &
                        grid(g) % n_nodes,  &
                        grid(g) % n_cells)
    call Save_Vtu_Faces(grid(g))

    ! Create 1D file (used for channel or pipe flow)
    call Probe_1d_Nodes(grid(g))

  end do

  end program
