!==============================================================================!
  subroutine Load_Cgns(grid)
!------------------------------------------------------------------------------!
!   Reads the Cgns mesh file.                                                  !
!   This function retrieves data as Load_Neu does, but for .cgns               !
!   https://cgns.github.io/index.html                                          !
!   Recomendations for Salome mesher:                                          !
!   use "group" option during export to .cgns                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod,  only: problem_name
  use Grid_Mod,  only: Grid_Type,                    &
                       Grid_Mod_Print_Bnd_Cond_List
  use Cgns_Mod   ! plenty of functions are used
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in
  integer           :: c, i, j, bc, base, block, sect, int, coord, mode
  integer           :: cgns_1, cgns_2, cgns_3, cgns_4, cgns_5, cell_type
  integer           :: parent_flag =1
!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = '.cgns'

  file_name = name_in

  !---------------------------------------!
  !                                       !
  !                                       !
  !   First run: just read info from DB   !
  !                                       !
  !                                       !
  !---------------------------------------!
  call Cgns_Mod_Initialize_Counters

  ! Open a CGNS file (->file_id)
  mode = CG_MODE_READ
  call Cgns_Mod_Open_File(file_name, mode)

  ! Read number of CGNS bases in file_id (->n_bases)
  call Cgns_Mod_Read_Number_Of_Bases_In_File

  !------------------------------!
  !                              !
  !   Browse through all bases   !
  !                              !
  !------------------------------!
  do base = 1, n_bases

    ! Read CGNS base information in base_id
    call Cgns_Mod_Read_Base_Info(base)

    !-------------------------------------------!
    !   Browse through all blocks in the base   !
    !-------------------------------------------!
    call Cgns_Mod_Read_Number_Of_Blocks_In_Base(base)

    do block = 1, cgns_base(base) % n_blocks

      ! Read block information
      ! Gives: n_nodes, n_cells
      call Cgns_Mod_Read_Block_Info(base, block)

      ! Read type of block
      ! Tells if structured or unstructured
      call Cgns_Mod_Read_Block_Type(base, block)

      ! Read number of interfacafical domains in block
      call Cgns_Mod_Read_Number_Of_Block_Interfaces(base, block)

      do int = 1, cgns_base(base) % block(block) % n_interfaces
        ! Reads info on interfaces
        call Cgns_Mod_Read_Interface_Info(base, block, int)
      end do ! interfaces

      ! Browse through all boundary conditions
      call Cgns_Mod_Read_Number_Of_Bnd_Conds_In_Block(base, block)

      do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
        ! Gives b.c. names and colours
        call Cgns_Mod_Read_Bnd_Conds_Info(base, block, bc)
      end do

      !-----------------------------------------!
      !   Browse through all element sections   !
      !-----------------------------------------!
      call Cgns_Mod_Read_Number_Of_Element_Sections(base, block)

      do sect = 1, cgns_base(base) % block(block) % n_sects

        ! Read info for an element section (including b.c.)
        ! Gives: cell_type, first_cell, last_cell
        call Cgns_Mod_Read_Section_Info(base, block, sect)

      end do ! elements sections

    end do ! blocks
  end do ! bases

  !------------------------!
  !                        !
  !   First run: results   !
  !                        !
  !------------------------!
  if(verbose) then
    print "(a)",     " #======================================================="
    print "(a)",     " # First run finished!"
    print "(a,i35)", " # - number of nodes: ",                      cnt_nodes
    print "(a,i35)", " # - number of cells: ",                      cnt_cells
    print "(a,i31)", " # - number of hex cells: ",                  cnt_hex
    print "(a,i26)", " # - number of pyramids cells: ",             cnt_pyr
    print "(a,i29)", " # - number of prism cells: ",                cnt_wed
    print "(a,i29)", " # - number of tetra cells: ",                cnt_tet
    print "(a,i29)", " # - number of mixed cells: ",                cnt_mix
    print "(a,i13)", " # - number of triangles faces on boundary: ",cnt_bnd_tri
    print "(a,i17)", " # - number of quads faces on boundary: ",    cnt_bnd_qua
    print "(a,i15)", " # - number of boundary conditions faces: ",  cnt_bnd_tri&
                                                                  + cnt_bnd_qua
    print "(a,i12)", " # - number of triangles faces on interface: ", &
                                                                    cnt_int_tri
    print "(a,i16)", " # - number of quads faces on interface: ", &
                                                                    cnt_int_qua
    print "(a,i21)", " # - number of boundary conditions: ",        cnt_bnd_cond
    if (cnt_bnd_tri + cnt_bnd_qua .eq. 0) then
      print "(a)",   " # No boundary faces were found !"
      stop
    end if
    print "(a)",    " #-------------------------------------------------------"
  end if

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  grid % n_nodes     = cnt_nodes
  grid % n_cells     = cnt_cells
  grid % n_bnd_cells = cnt_bnd_tri + cnt_bnd_qua

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  grid % n_bnd_cond  = cnt_bnd_cond
  allocate(grid % bnd_cond % name(cnt_bnd_cond))
  do i = 1, cnt_bnd_cond
    call To_Upper_Case( bnd_cond_names(i) )
    grid % bnd_cond % name(i) = bnd_cond_names(i)
  end do

  !----------------!
  !   Interfaces   !
  !----------------!
  ! Explained in Cgns_Mod_Merge_Nodes
  allocate(interface_cells(2, cnt_int_qua + cnt_int_tri, 4, cnt_int))
  interface_cells = -1

  call Allocate_Memory(grid)

  !-------------------------------------!
  !                                     !
  !                                     !
  !   Second run: read arrays from DB   !
  !                                     !
  !                                     !
  !-------------------------------------!
  call Cgns_Mod_Initialize_Counters

  print *, '# Filling arrays..'

  !------------------------------!
  !                              !
  !   Browse through all bases   !
  !                              !
  !------------------------------!
  do base = 1, n_bases

    !-------------------------------!
    !   Browse through all blocks   !
    !-------------------------------!
    do block = 1, cgns_base(base) % n_blocks
      print "(a,a21)", " # Block name: ",  &
               trim(cgns_base(base) % block(block) % name)

      ! Count block, just for information
      cnt_blocks = cnt_blocks + 1

      !----------------------------!
      !   Read coordinates block   !
      !----------------------------!

      ! Reads number of coordinates arrays from block
      ! Essentially, reads number 3
      call Cgns_Mod_Read_Number_Of_Coordinates_In_Block(base, block)

      ! Read x, y and z coordinates
      do coord = 1, cgns_base(base) % block (block) % n_coords

        ! Reads coord_name of coord_id(->coord_name)
        call Cgns_Mod_Read_Coordinate_Info(base, block, coord)

        ! Read grid coordinates (-> x_, y_, z_coord)
        call Cgns_Mod_Read_Coordinate_Array(base, block, coord, grid)

      end do ! coordinates

      !----------------------!
      !   Read cells block   !
      !----------------------!

      ! Browse through all sections to read elements
      do sect = 1, cgns_base(base) % block(block) % n_sects

        ! Read element data (HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3/MIXED)
        call Cgns_Mod_Read_Section_Connections(base, block, sect, grid)

      end do ! elements sections

      cnt_nodes = cnt_nodes + cgns_base(base) % block(block) % mesh_info(1)
      cnt_cells = cnt_cells + cgns_base(base) % block(block) % mesh_info(2)

    end do ! blocks
  end do ! bases

  !------------------!
  !   Find parents   !
  !------------------!
  !if(parent_flag .eq. 0) then
  !  call Find_Parents(grid)
  !end if

  print "(a)",     " #-------------------------------------------------"
  print "(a,i13)", " # Total number of nodes:             ", cnt_nodes
  print "(a,i13)", " # Total number of cells:             ", cnt_cells
  print "(a,i13)", " # Total number of blocks:            ", cnt_blocks
  print "(a,i13)", " # Total number of boundary sections: ", grid % n_bnd_cond
  print "(a,i13)", " # Total number of boundary cells:    ", grid % n_bnd_cells
  print "(a)",     " #-------------------------------------------------"

  !---------------------!
  !   Merge the nodes   !
  !---------------------!
  if(cnt_blocks .gt. 1) then
    call Cgns_Mod_Merge_Nodes(grid)
  end if

  !---------------------------------!
  !   Read block (material?) data   !
  !---------------------------------!
  grid % material % name = "AIR"

  !-----------------------------------------------------------------!
  !   Correct boundary conditions directions for hexahedral cells   !
  !   (They are not the same in CGNS and Gambit's neutral format.)  !
  !-----------------------------------------------------------------!
  if(verbose) then
    print "(a)", " # HEX_8 colors (sample):"
    print "(a, 7a7)", " # ", &
      "c", "bnd_c1", "bnd_c2", "bnd_c3", "bnd_c4", "bnd_c5", "bnd_c6"
    i = 0
  end if

  !if(parent_flag .ne. 0) then
    cnt_bnd_cells = 0
    do base = 1, n_bases
      do block = 1, cgns_base(base) % n_blocks
        do sect = 1, cgns_base(base) % block(block) % n_sects
          cell_type = cgns_base(base) % block(block) % section(sect) % cell_type

          if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) then

            do c = cgns_base(base) % block(block) % section(sect) % first_cell,&
                   cgns_base(base) % block(block) % section(sect) % last_cell
              cgns_1 = grid % cells_bnd_color(1,c)
              cgns_2 = grid % cells_bnd_color(2,c)
              cgns_3 = grid % cells_bnd_color(3,c)
              cgns_4 = grid % cells_bnd_color(4,c)
              cgns_5 = grid % cells_bnd_color(5,c)
              grid % cells_bnd_color(4,c) = cgns_5
              grid % cells_bnd_color(3,c) = cgns_4
              grid % cells_bnd_color(2,c) = cgns_3
              grid % cells_bnd_color(1,c) = cgns_2
              grid % cells_bnd_color(5,c) = cgns_1

              do j = 1, 6
                if( grid % cells_bnd_color(j,c) .ne. 0 ) then
                  cnt_bnd_cells = cnt_bnd_cells + 1
                end if
              end do

              if(verbose .and. i < 7) then
                print "(a,7i7)"," # ",c, (grid % cells_bnd_color(j,c), j = 1, 6)
                i = i + 1
              end if

            end do ! c
          end if ! ElementTypeName(cell_type) .eq. 'HEXA_8'

        end do ! elements sections
      end do ! blocks
    end do ! bases
    print "(a)", " #--------------------------------------------------"
    print "(a,i19)", " # Corrected hex boundary cells: ", cnt_bnd_cells
  !end if

  call Grid_Mod_Print_Bnd_Cond_List(grid)

  end subroutine
