!==============================================================================!
  subroutine Load_Gambit(grid, file_name)
!------------------------------------------------------------------------------!
!   Reads the Gambits neutral file format.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  character(SL)   :: file_name
!-----------------------------------[Locals]-----------------------------------!
  integer             :: i, j, n_blocks, n_bnd_sect, dum1, dum2, fu
  integer,allocatable :: temp(:)
  integer             :: c, dir
!==============================================================================!

  call File_Mod_Open_File_For_Reading(file_name, fu)

  !------------------------------------------!
  !   Gambit can't handle polyhedral grids   !
  !------------------------------------------!
  grid % polyhedral = .false.

  !------------------------!
  !   Browse over header   !
  !------------------------!
  do i = 1, 10
    call File_Mod_Read_Line(fu)

    ! At line containing: "NUMNP NELEM NGRPS NBSETS NDFCD NDFVL" jump out
    if(line % n_tokens .eq. 6) then
      if(line % tokens(1) .eq. 'NUMNP') then
        goto 1
      end if
    end if
  end do

  !--------------------------------------------------!
  !   Read the first line with usefull information   !
  !--------------------------------------------------!
1 continue
  call File_Mod_Read_Line(fu)

  read(line % tokens(1),*) grid % n_nodes
  read(line % tokens(2),*) grid % n_cells
  read(line % tokens(3),*) n_blocks
  read(line % tokens(4),*) n_bnd_sect

  print '(a38,i9)', '# Total number of nodes:             ', grid % n_nodes
  print '(a38,i9)', '# Total number of cells:             ', grid % n_cells
  print '(a38,i9)', '# Total number of blocks:            ', n_blocks
  print '(a38,i9)', '# Total number of boundary sections: ', n_bnd_sect

  !------------------------------!
  !   Count the boundary cells   !
  !------------------------------!
  grid % n_bnd_cells = 0
  do
    call File_Mod_Read_Line(fu)
    if( line % tokens(1) .eq. 'BOUNDARY' ) then
      do j = 1, n_bnd_sect
        if(j>1) call File_Mod_Read_Line(fu) ! BOUNDARY CONDITIONS
        call File_Mod_Read_Line(fu)
        read(line % tokens(3),*) dum1
        grid % n_bnd_cells = grid % n_bnd_cells + dum1 
        do i = 1, dum1
          read(fu,*) c, dum2, dir
        end do
        call File_Mod_Read_Line(fu)         ! ENDOFSECTION
      end do
      print '(a38,i9)', '# Total number of boundary cells:    ',  &
            grid % n_bnd_cells
      goto 2
    end if
  end do

2 rewind(fu)

  !------------------------!
  !   Browse over header   !
  !------------------------!
  do i = 1, 10
    call File_Mod_Read_Line(fu)

    ! At line containing: "ENDOFSECTION" jump out
    if(line % n_tokens .eq. 1) then
      if(line % tokens(1) .eq. 'ENDOFSECTION') then
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
  call Allocate_Memory(grid)

  allocate(temp(grid % n_cells)); temp=0

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  call File_Mod_Read_Line(fu)          ! NODAL COORDINATES
  do i = 1, grid % n_nodes
    call File_Mod_Read_Line(fu)
    read(line % tokens(2),*) grid % xn(i)
    read(line % tokens(3),*) grid % yn(i)
    read(line % tokens(4),*) grid % zn(i)
  end do
  call File_Mod_Read_Line(fu)          ! ENDOFSECTION

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  call File_Mod_Read_Line(fu)          ! ELEMENTS/CELLS
  do i = 1, grid % n_cells
    read(fu,'(i8,1x,i2,1x,i2,1x,7i8:/(15x,7i8:))') dum1, dum2, &
           grid % cells_n_nodes(i),                           &
          (grid % cells_n(j,i), j = 1, grid % cells_n_nodes(i))

    ! Nodes 3,4 and 7,8 should be swapped for hexahedral
    ! cells to be compatible with writing in vtu format
    if( grid % cells_n_nodes(i) .eq. 8 ) then
      call Swap_Int(grid % cells_n(3, i), grid % cells_n(4, i))
      call Swap_Int(grid % cells_n(7, i), grid % cells_n(8, i))
    end if

    ! Nodes 3, 4 should be swapped for pyramids
    ! to be compatible with files in vtu format
    if( grid % cells_n_nodes(i) .eq. 5 ) then
      call Swap_Int(grid % cells_n(3, i), grid % cells_n(4, i))
    end if

  end do
  call File_Mod_Read_Line(fu)          ! ENDOFSECTION

  !----------------------!
  !   Material section   !
  !----------------------!
  do j = 1, n_blocks
    call File_Mod_Read_Line(fu)               ! ELEMENT GROUP
    call File_Mod_Read_Line(fu)
    read(line % tokens(4),'(i10)') dum1       ! number of cells in this group
    call File_Mod_Read_Line(fu)               ! block*
    call File_Mod_Read_Line(fu)               ! 0
    read(fu,'(10i8)') (temp(i), i = 1, dum1)  ! read all cells in the group
    call File_Mod_Read_Line(fu)               ! ENDOFSECTION
  end do

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  grid % n_bnd_cond = n_bnd_sect
  allocate(grid % bnd_cond % name(n_bnd_sect))

  do j = 1, n_bnd_sect
    call File_Mod_Read_Line(fu)        ! BOUNDARY CONDITIONS
    call File_Mod_Read_Line(fu)
    call To_Upper_Case(  line % tokens(1)  )
    grid % bnd_cond % name(j) = line % tokens(1)
    read(line % tokens(3),'(i8)') dum1
    do i = 1, dum1
      read(fu,*) c, dum2, dir
      grid % cells_bnd_color(dir,c) = j
    end do
    call File_Mod_Read_Line(fu)        ! ENDOFSECTION
  end do

  !------------------------------------!
  !   Pring boundary conditions info   !
  !------------------------------------!
  call Grid_Mod_Print_Bnd_Cond_List(grid)

  close(fu)

  end subroutine
