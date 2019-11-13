!==============================================================================!
  subroutine Save_Cgns_Cells(grid, sub)
!------------------------------------------------------------------------------!
!   Writes in 3-D unstructured grid to files 'file_name.cgns'                  !
!   Valid for both parallel and seqential access                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  integer           :: sub
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out
  integer           :: c, base, block, sect, coord
!==============================================================================!

  !-----------------------!
  !                       !
  !   Create .cgns file   !
  !                       !
  !-----------------------!

  if (.not. mesh_written) then
    ! problem_name.cgns
    call File_Mod_Set_Name(name_out, extension='.cgns')
    if (sub .lt. 2) print *, '# Creating the file with mesh  : ', trim(name_out)
  else
    ! problem_name-ts??????.cgns
    name_out = trim(file_name)
    if (sub .lt. 2) print *, '# Creating the file with fields: ', trim(name_out)
  end if

  file_mode = CG_MODE_WRITE
  call Cgns_Mod_Open_File(name_out, file_mode)

  call Cgns_Mod_Initialize_Counters

  ! Count number of 3d cell type elements
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % cells_n_nodes(c) .eq. 8) cnt_hex = cnt_hex + 1
    if(grid % cells_n_nodes(c) .eq. 6) cnt_wed = cnt_wed + 1
    if(grid % cells_n_nodes(c) .eq. 5) cnt_pyr = cnt_pyr + 1
    if(grid % cells_n_nodes(c) .eq. 4) cnt_tet = cnt_tet + 1
  end do

  !-----------------!
  !                 !
  !   Bases block   !
  !                 !
  !-----------------!
  n_bases = 1
  allocate(cgns_base(n_bases))

  base = 1
  cgns_base(base) % name = "Base 1"
  cgns_base(base) % cell_dim = 3
  cgns_base(base) % phys_dim = 3

  call Cgns_Mod_Write_Base_Info(base)

  !-----------------!
  !                 !
  !   Zones block   !
  !                 !
  !-----------------!

  cgns_base(base) % n_blocks = 1
  allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

  block = 1
  cgns_base(base) % block(block) % name = "Zone 1"
  c = grid % n_nodes
  call Comm_Mod_Global_Sum_Int(c)
  cgns_base(base) % block(block) % mesh_info(1) = c
  c = grid % n_cells - grid % comm % n_buff_cells
  call Comm_Mod_Global_Sum_Int(c)
  cgns_base(base) % block(block) % mesh_info(2) = c
  cgns_base(base) % block(block) % mesh_info(3) = 0

  call Cgns_Mod_Write_Block_Info(base, block)


  !-----------------------!
  !                       !
  !   Coordinates block   !
  !                       !
  !-----------------------!

  coord = 1
  cgns_base(base) % block(block) % coord_name(coord) = "CoordinateX"

  coord = 2
  cgns_base(base) % block(block) % coord_name(coord) = "CoordinateY"

  coord = 3
  cgns_base(base) % block(block) % coord_name(coord) = "CoordinateZ"

  ! Actually write grid coordinates in DB
  if (.not. mesh_written) then
    do coord = 1, cgns_base(base) % cell_dim
      call Cgns_Mod_Write_Coordinate_Array(base, block, coord, grid)
    end do
  end if

  !-----------------------------!
  !                             !
  !   Cells connections block   !
  !                             !
  !-----------------------------!

  cgns_base(base) % block(block) % n_sects = 4
  allocate(cgns_base(base) % block(block) % section( &
    cgns_base(base) % block(block) % n_sects))

  sect = 1
  cgns_base(base) % block(block) % section(sect) % name = "Hexagons"
  cgns_base(base) % block(block) % section(sect) % cell_type = HEXA_8
  cgns_base(base) % block(block) % section(sect) % first_cell = 1
  cgns_base(base) % block(block) % section(sect) % last_cell = cnt_hex

  sect = 2
  cgns_base(base) % block(block) % section(sect) % name = "Wedges"
  cgns_base(base) % block(block) % section(sect) % cell_type = PENTA_6
  cgns_base(base) % block(block) % section(sect) % first_cell = 1
  cgns_base(base) % block(block) % section(sect) % last_cell = cnt_wed

  sect = 3
  cgns_base(base) % block(block) % section(sect) % name = "Pyramids"
  cgns_base(base) % block(block) % section(sect) % cell_type = PYRA_5
  cgns_base(base) % block(block) % section(sect) % first_cell = 1
  cgns_base(base) % block(block) % section(sect) % last_cell = cnt_pyr

  sect = 4
  cgns_base(base) % block(block) % section(sect) % name = "Tetrahedrons"
  cgns_base(base) % block(block) % section(sect) % cell_type = TETRA_4
  cgns_base(base) % block(block) % section(sect) % first_cell = 1
  cgns_base(base) % block(block) % section(sect) % last_cell = cnt_tet
  
  ! actually write grid connections in DB
  if (.not. mesh_written) then

    do sect = 1, cgns_base(base) % block(block) % n_sects
      call Cgns_Mod_Write_Section_Connections(base, block, sect, grid)
    end do

  else ! write an link to actual mesh inside a different DB file
    call Write_Link_To_Mesh_In_File(file_with_mesh, base, block)
  end if

  ! Close DB
  call Cgns_Mod_Close_File

  if (.not. mesh_written .and. sub .lt. 2) &
    print *, '# Have written unstructured grid to file: ',trim(name_out)

  deallocate(cgns_base)

  end subroutine
