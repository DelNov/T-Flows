#include "../Shared/Assert.h90"

!==============================================================================!
  program Convert_Prog
!------------------------------------------------------------------------------!
!>  Convert's main function.  Its primary function is to import grid files in
!>  various formats, perform necessary transformations and calculations, and
!>  prepare them for CFD simulations with T-Flows. It handles various file
!>  formats such as FLUENT, GAMBIT, GMSH, and OBJ.
!------------------------------------------------------------------------------!
!   Program Flow                                                               !
!                                                                              !
!   * Initialization: Begins with profiling, followed by the display of the    !
!     program's logo.                                                          !
!   * File name and format: The user is prompted to enter the grid file name.  !
!     The program then guesses the file's format for appropriate processing.   !
!   * Grid processing: Depending on the detected file format (FLUENT, GAMBIT,  !
!     GMSH, OBJ), the corresponding loading procedure is invoked from          !
!     Convert_Mod.                                                             !
!   * Cityscape handling: If the grid represents a cityscape (determined from  !
!     the file name), specific city-related processing is performed.           !
!   * Dual grid option: The user is given the choice to create a dual grid.    !
!     If selected, both primal and dual grids are processed.                   !
!   * Geometrical calculations and sorting: For each grid, the program         !
!     calculates geometrical properties, sorts cells and faces, and performs   !
!     other preparatory operations.                                            !
!   * Visualization and post-processing data saving: The grids are saved in    !
!     various formats (.cfn and .dim for further processing with T-Flows       !
!     and .vtu for visualization).                                             !
!   * 1D probe creation: For certain grid types, 1D probe files are generated  !
!     for specific flow analyses.                                              !
!   * Final Steps: The program prints grid statistics and concludes with       !
!      profiling data finalization.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Convert_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include '../Shared/Probe_1d_Nodes.h90'
  end interface
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG    = .false.
  logical, parameter :: SAVE_CAS = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type)     :: Grid(2)        ! grid to be converted and its dual
  type(Porosity_Type) :: Por            ! porosity
  character(SL)       :: answer
  character(SL)       :: app_up
  character(SL)       :: file_name
  character(SL)       :: file_format    ! 'UNKNOWN', 'FLUENT', 'GAMBIT', 'GMSH'
  integer             :: s, l, p, g, n_grids
  logical             :: city, porosity
!==============================================================================!

  ! Initialize program profler
  call Profiler % Start('Main')

  ! Open with a logo
  call Convert % Logo_Con()

  print '(a)', ' #========================================================'
  print '(a)', ' # Enter the grid file name you are importing (with ext.):'
  print '(a)', ' #--------------------------------------------------------'
  file_name = File % Single_Word_From_Keyboard()

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
  if(l .ge. 10) then
    app_up = problem_name(1)(l-7:l-4)
    call String % To_Upper_Case(app_up)
    if(app_up .eq. 'CITY') city = .true.
  end if

  porosity = .false.
  if(l .ge. 15) then
    app_up = problem_name(1)(l-11:l-4)
    call String % To_Upper_Case(app_up)
    if(app_up .eq. 'POROSITY') porosity = .true.
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
    call Grid(1) % Save_Vtu_Faces((/0, 0/))

    print '(a)', ' #========================================================'
    print '(a)', ' # Enter STL forrest file name to plant trees (with ext.):'
    print '(a)', ' #--------------------------------------------------------'
    read(*,*) file_name

    ! Read the forrest STL file and plant the trees
    call Convert % Load_Forrest(Grid, file_name)
    call Grid(2) % Save_Vtu_Faces((/0, 0/))

    ! Finalize program profler
    call Profiler % Stop('Main')
    call Profiler % Statistics(indent=1)
    stop
  end if

  if(city .or. porosity) then
    call Convert % Find_Boundary_Faces(Grid(1))
    call Convert % Find_Inside_Faces  (Grid(1))
  end if

  ! Take of the city-related domains
  if(city) then
    call Convert % Insert_Buildings(Grid(1))
  end if

  if(porosity) then
    call Grid(1) % Calculate_Cell_Centers()
    call Work % Allocate_Work(Grid)
    call Por % Create_Porosity(Grid(1))
    call Work % Finalize_Work()
    call Convert % Remove_Porosities(Grid(1))
  end if

  ! For Gambit and Gmsh grids, no face information is stored
  if(file_format .eq. 'GAMBIT' .or. file_format .eq. 'GMSH') then
    call Convert % Find_Boundary_Faces(Grid(1))
    call Convert % Find_Inside_Faces  (Grid(1))
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
  answer = File % Single_Word_From_Keyboard()
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
    call Convert % Calculate_Geometry(Grid(g),         &
                                      (g.eq.n_grids),  &  ! should you scale
                                      g)

    ! Keep in mind that Grid_Mod_Calculate_Wall_Distance is ...
    ! ... faster if it is called after Grid_Mod_Sort_Faces_By_Region
    call Grid(g) % Sort_Cells_By_Coordinates()

    call Grid(g) % Sort_Faces_By_Region()
    if( (g-n_grids) .eq. 0) then
      call Grid(g) % Calculate_Wall_Distance()
    end if
    call Grid(g) % Find_Cells_Faces()
    if(DEBUG) call Grid(g) % Check_Cells_Closure()

    call Grid(g) % Initialize_New_Numbers()

    ! Note #1 about shadows:
    ! At this point you have Grid % n_faces faces and Grid % n_shadows (on top)
    ! and they are pointing to each other.  Besides, both real face and its
    ! shadow have the same c1 and c2, both inside cells with positive indices
    ! Real faces which do not have shadows have index "0" for shadow.
    do s = 1, Grid(g) % n_faces + Grid(g) % n_shadows
      if(Grid(g) % faces_s(s) .ne. 0) then
        Assert(Grid(g) % faces_s(Grid(g) % faces_s(s)) .eq. s)
      end if
    end do
    ! Similar note is in Generate, also called Note #1

    call Grid(g) % Print_Grid_Statistics()

    !-------------------------------!
    !   Save files for processing   !
    !-------------------------------!
    call Grid(g) % Save_Cfn((/0, 0/),               &
                            Grid(g) % n_nodes,      &
                            Grid(g) % n_cells,      &
                            Grid(g) % n_faces,      &
                            Grid(g) % n_shadows,    &
                            Grid(g) % n_bnd_cells)

    call Grid(g) % Save_Dim((/0, 0/))

    !-----------------------------------------------------!
    !   Save grid for visualisation and post-processing   !
    !-----------------------------------------------------!

    ! Create output in vtu format
    call Grid(g) % Save_Vtu_Cells((/0, 0/),           &
                                  Grid(g) % n_nodes,  &
                                  Grid(g) % n_cells)
    call Grid(g) % Save_Vtu_Faces((/0, 0/))
    call Grid(g) % Save_Vtu_Faces((/0, 0/), plot_inside=.false.)
    call Grid(g) % Save_Vtu_Faces((/0, 0/), plot_shadows=.true.)

    ! Save Fluent's .cas file
    if(SAVE_CAS) then
      call Convert % Save_Fluent(Grid(g))
    end if

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
