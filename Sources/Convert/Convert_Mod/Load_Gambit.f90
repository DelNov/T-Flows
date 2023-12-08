!==============================================================================!
  subroutine Load_Gambit(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!>  Subroutine is designed to read grid files in Gambit's neutral file format
!>  and populate the provided Grid object with the necessary data.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File reading: Opens the specified Gambit file in ASCII mode for reading. !
!   * Header processing: Scans through the file's header to extract key        !
!     information like the number of nodes, cells, blocks, and bnd sections.   !
!   * Boundary cells counting: Counts the number of boundary cells by reading  !
!     relevant sections in the file.                                           !
!   * Memory allocation: Allocates memory for various arrays within the Grid   !
!     object based on the extracted information.                               !
!   * Reading nodal coordinates: Reads coordinates for each node from the file !
!     and stores them in the Grid object.                                      !
!   * Reading cell nodes: Reads and stores the nodes of each cell. Special     !
!     handling is done for hexahedral cells and pyramids to ensure             !
!     compatibility with .vtu format conventions (swapping certain nodes).     !
!   * Material section processing: Handles the material section in the Gambit  !
!     file. Each block’s material information is read, but the details of the  !
!     processing are not explicitly stored.                                    !
!   * Boundary conditions processing: Reads boundary condition sections,       !
!     including the name of each boundary condition and the cells associated   !
!     with each condition. The information is stored in the Grid object.       !
!   * Finalization: Closes the file and stops profiling the subroutine.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid       !! grid being converted
  character(SL)       :: file_name  !! file name
!-----------------------------------[Locals]-----------------------------------!
  integer             :: i, j, n_blocks, n_bnd_sect, dum1, dum2, fu
  integer,allocatable :: temp(:)
  integer             :: c, dir, nodes(8)
!==============================================================================!

  call Profiler % Start('Load_Gambit')

  call File % Open_For_Reading_Ascii(file_name, fu)

  !------------------------------------------!
  !   Gambit can't handle polyhedral grids   !
  !------------------------------------------!
  Grid % polyhedral = .false.

  !------------------------!
  !   Browse over header   !
  !------------------------!
  do i = 1, 10
    call File % Read_Line(fu)

    ! At line containing: "NUMNP NELEM NGRPS NBSETS NDFCD NDFVL" jump out
    if(Line % n_tokens .eq. 6) then
      if(Line % tokens(1) .eq. 'NUMNP') then
        goto 1
      end if
    end if
  end do

  !--------------------------------------------------!
  !   Read the first line with usefull information   !
  !--------------------------------------------------!
1 continue
  call File % Read_Line(fu)

  read(Line % tokens(1),*) Grid % n_nodes
  read(Line % tokens(2),*) Grid % n_cells
  read(Line % tokens(3),*) n_blocks
  read(Line % tokens(4),*) n_bnd_sect

  print '(a38,i9)', '# Total number of nodes:             ', Grid % n_nodes
  print '(a38,i9)', '# Total number of cells:             ', Grid % n_cells
  print '(a38,i9)', '# Total number of blocks:            ', n_blocks
  print '(a38,i9)', '# Total number of boundary sections: ', n_bnd_sect

  !------------------------------!
  !   Count the boundary cells   !
  !------------------------------!
  Grid % n_bnd_cells = 0
  do
    call File % Read_Line(fu)
    if( Line % tokens(1) .eq. 'BOUNDARY' ) then
      do j = 1, n_bnd_sect
        if(j>1) call File % Read_Line(fu) ! BOUNDARY CONDITIONS
        call File % Read_Line(fu)
        read(Line % tokens(3),*) dum1
        Grid % n_bnd_cells = Grid % n_bnd_cells + dum1 
        do i = 1, dum1
          read(fu,*) c, dum2, dir
        end do
        call File % Read_Line(fu)         ! ENDOFSECTION
      end do
      print '(a38,i9)', '# Total number of boundary cells:    ',  &
            Grid % n_bnd_cells
      goto 2
    end if
  end do

2 rewind(fu)

  !------------------------!
  !   Browse over header   !
  !------------------------!
  do i = 1, 10
    call File % Read_Line(fu)

    ! At line containing: "ENDOFSECTION" jump out
    if(Line % n_tokens .eq. 1) then
      if(Line % tokens(1) .eq. 'ENDOFSECTION') then
        goto 3
      end if
    end if
  end do

3 continue

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  call Convert % Allocate_Memory(Grid)

  allocate(temp(Grid % n_cells)); temp=0

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  call File % Read_Line(fu)          ! NODAL COORDINATES
  do i = 1, Grid % n_nodes
    call File % Read_Line(fu)
    read(Line % tokens(2),*) Grid % xn(i)
    read(Line % tokens(3),*) Grid % yn(i)
    read(Line % tokens(4),*) Grid % zn(i)
  end do
  call File % Read_Line(fu)          ! ENDOFSECTION

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  call File % Read_Line(fu)          ! ELEMENTS/CELLS
  do i = 1, Grid % n_cells
    read(fu,'(i8,1x,i2,1x,i2,1x,7i8:/(15x,7i8:))') dum1, dum2, &
           Grid % cells_n_nodes(i),                           &
          (nodes(j), j = 1, Grid % cells_n_nodes(i))

    call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,Grid % cells_n_nodes(i)/))
    do j = 1, Grid % cells_n_nodes(i)
      Grid % cells_n(j, i) = nodes(j)
    end do

    ! Nodes 3,4 and 7,8 should be swapped for hexahedral
    ! cells to be compatible with writing in vtu format
    if( Grid % cells_n_nodes(i) .eq. 8 ) then
      call Swap_Int(Grid % cells_n(3, i), Grid % cells_n(4, i))
      call Swap_Int(Grid % cells_n(7, i), Grid % cells_n(8, i))
    end if

    ! Nodes 3, 4 should be swapped for pyramids
    ! to be compatible with files in vtu format
    if( Grid % cells_n_nodes(i) .eq. 5 ) then
      call Swap_Int(Grid % cells_n(3, i), Grid % cells_n(4, i))
    end if

  end do
  call File % Read_Line(fu)          ! ENDOFSECTION

  !----------------------!
  !   Material section   !
  !----------------------!
  do j = 1, n_blocks
    call File % Read_Line(fu)                 ! ELEMENT GROUP
    call File % Read_Line(fu)
    read(Line % tokens(4),'(i10)') dum1       ! number of cells in this group
    call File % Read_Line(fu)                 ! block*
    call File % Read_Line(fu)                 ! 0
    read(fu,'(10i8)') (temp(i), i = 1, dum1)  ! read all cells in the group
    call File % Read_Line(fu)                 ! ENDOFSECTION
  end do

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  call Grid % Allocate_Regions(n_bnd_sect)

  do j = 1, n_bnd_sect
    call File % Read_Line(fu)        ! BOUNDARY CONDITIONS
    call File % Read_Line(fu)
    call String % To_Upper_Case(Line % tokens(1))
    Grid % region % name(j) = Line % tokens(1)
    read(Line % tokens(3),'(i8)') dum1
    do i = 1, dum1
      read(fu,*) c, dum2, dir
      Grid % cells_bnd_region(dir,c) = j
    end do
    call File % Read_Line(fu)        ! ENDOFSECTION
  end do

  !------------------------------------!
  !   Pring boundary conditions info   !
  !------------------------------------!
  call Grid % Print_Regions_List()

  close(fu)

  call Profiler % Stop('Load_Gambit')

  end subroutine
