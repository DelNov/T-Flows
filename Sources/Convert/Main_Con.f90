#include "../Shared/Assert.h90"

!==============================================================================!
  program Convert_Prog
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Convert_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include '../Shared/Probe_1d_Nodes.h90'
  end interface
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: Grid(2)        ! grid to be converted and its dual
  character(SL)   :: answer
  character(SL)   :: app_up
  character(SL)   :: file_name
  character(SL)   :: file_format    ! 'UNKNOWN', 'FLUENT', 'GAMBIT', 'GMSH'
  integer         :: l, p, g, n_grids
  logical         :: city
!==============================================================================!

  ! Initialize program profler
  call Profiler % Start('Main')

  ! Open with a logo
  call Convert % Logo_Con()

  print '(a)', ' #========================================================'
  print '(a)', ' # Enter the grid file name you are importing (with ext.):'
  print '(a)', ' #--------------------------------------------------------'
  read(*,*) file_name

  !-----------------------------------------------!
  !                                               !
  !   Make an educated guess of the file format   !
  !                                               !
  !-----------------------------------------------!
  call Convert % Guess_Format(file_name, file_format)

  !-------------------------------------------!
  !                                           !
  !   Extract the problem name and store it   !
  !                                           !
  !-------------------------------------------!
  l = len_trim(file_name)
  p = index(file_name(1:l), '.', back=.true.)

  problem_name(1) = file_name(1:p-1)

  Grid(1) % name = problem_name(1)
  call String % To_Upper_Case(Grid(1) % name)
  city = .false.
  if(l > 10) then
    app_up = problem_name(1)(l-7:l-4)
    call String % To_Upper_Case(app_up)
    if(app_up .eq. 'CITY') city = .true.
  end if

  !----------------------------------------!
  !                                        !
  !   Read the file and start conversion   !
  !                                        !
  !----------------------------------------!
  if(file_format .eq. 'FLUENT') then
    call Convert % Load_Fluent(Grid(1), file_name)
  end if
  if(file_format .eq. 'GAMBIT') then
    call Convert % Load_Gambit(Grid(1), file_name)
  end if
  if(file_format .eq. 'GMSH') then
    call Convert % Load_Gmsh(Grid(1), file_name)
    call Convert % Find_Parents(Grid(1))
  end if
  if(file_format .eq. 'OBJ') then

    ! Read the single tree
    call Convert % Load_Obj(Grid(1), file_name)
    call Grid(1) % Save_Vtu_Faces()

    print '(a)', ' #========================================================'
    print '(a)', ' # Enter STL forrest file name to plant trees (with ext.):'
    print '(a)', ' #--------------------------------------------------------'
    read(*,*) file_name

    ! Read the forrest STL file and plant the trees
    call Convert % Load_Forrest(Grid, file_name)
    call Grid(2) % Save_Vtu_Faces()

    ! Finalize program profler
    call Profiler % Stop('Main')
    call Profiler % Statistics(indent=1)
    stop
  end if

  ! Sort cells in height first thing after reading	    
  if(city) then
    call Convert % Insert_Buildings(Grid(1))
  end if

  ! For Gambit and Gmsh grids, no face information is stored
  if(file_format .eq. 'GAMBIT' .or. file_format .eq. 'GMSH') then
    call Convert % Grid_Topology(Grid(1))
    call Convert % Find_Faces   (Grid(1))
  end if

  ! Some mesh generators (Gmsh for sure) can leave duplicate
  ! nodes in the grid. Check it and eliminate them with this
  if(file_format .eq. 'GMSH') then
    call Grid(1) % Merge_Duplicate_Nodes()
  end if

  !---------------------------------------------------!
  !                                                   !
  !   Decide if you are going for dual grid as well   !
  !                                                   !
  !---------------------------------------------------!
  print '(a)', ' #================================================='
  print '(a)', ' # Would you like to create a dual grid? (yes/no)'
  print '(a)', ' #-------------------------------------------------'
  call File % Read_Line(5)
  answer = Line % tokens(1)
  call String % To_Upper_Case(answer)

  n_grids = 1
  if(answer .eq. 'YES') n_grids = 2

  !-------------------------------!
  !                               !
  !   Browse through both grids   !
  !                               !
  !-------------------------------!
  do g = 1, n_grids

    if(n_grids .eq. 2) then
      print '(a)',              ' #======================================'
      print '(a)',              ' #                                      '
      if(g .eq. 1) print '(a)', ' # Processing the first (primal) grid'
      if(g .eq. 2) print '(a)', ' # Processing the second (dual) grid'
      print '(a)',              ' #                                    '
      print '(a)',              ' #--------------------------------------'
      if(g .eq. 2) call Convert % Create_Dual(Grid(1), Grid(2))
    end if

    !--------------------------------------!
    !   Calculate geometrical quantities   !
    !--------------------------------------!
    call Convert % Calculate_Geometry(Grid(g), g-n_grids)  ! if zero, ask

    ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
    ! ... faster if it is called after Grid_Mod_Sort_Faces_Smart
    call Grid(g) % Sort_Cells_Smart       ()
    call Grid(g) % Sort_Faces_Smart       ()
    if( (g-n_grids) .eq. 0) then
      call Grid(g) % Calculate_Wall_Distance()
    end if
    call Grid(g) % Find_Cells_Faces       ()

    call Grid(g) % Initialize_New_Numbers()

    ! Note #1 about shadows:
    ! At this point you have grid % n_faces faces and grid % n_shadows (on top)
    ! and they are pointing to each other.  Besides, both real face and its
    ! shadow have the same c1 and c2, both inside cells with positive indices
    ! Real faces which do not have shadows have index "0" for shadow.
    ! Should check like this: do s = 1, grid % n_faces + grid % n_shadows
    ! Should check like this:   write(20, '(99i9)') s, grid % faces_s(s)
    ! Should check like this: end do
    ! Similar note is in Generate, also called Note #1

    call Grid(g) % Print_Grid_Statistics()

    !-------------------------------!
    !   Save files for processing   !
    !-------------------------------!
    call Grid(g) % Save_Cfn(0,                      &
                            Grid(g) % n_nodes,      &
                            Grid(g) % n_cells,      &
                            Grid(g) % n_faces,      &
                            Grid(g) % n_shadows,    &
                            Grid(g) % n_bnd_cells)

    call Grid(g) % Save_Dim(0)

    !-----------------------------------------------------!
    !   Save grid for visualisation and post-processing   !
    !-----------------------------------------------------!

    ! Create output in vtu format
    call Grid(g) % Save_Vtu_Cells(0,                  &
                                  Grid(g) % n_nodes,  &
                                  Grid(g) % n_cells)
    call Grid(g) % Save_Vtu_Faces()
    call Grid(g) % Save_Vtu_Faces(plot_shadows=.true.)

    if( (g-n_grids) .eq. 0) then
      ! Create a template control file for this domain
      call Grid(g) % Write_Template_Control_File()

      ! Create 1D file (used for channel or pipe flow)
      call Probe_1d_Nodes(Grid(g))
    end if

  end do

  ! Finalize program profler
  call Profiler % Stop('Main')
  call Profiler % Statistics(indent=1)

  end program
