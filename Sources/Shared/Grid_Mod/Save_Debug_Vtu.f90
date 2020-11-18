!==============================================================================!
  subroutine Grid_Mod_Save_Debug_Vtu(grid, append,                           &
                                     scalar_cell, scalar_node, scalar_name,  &
                                     vector_cell, vector_node, vector_name)
!------------------------------------------------------------------------------!
!   Writes: name.vtu, name.faces.vtu, name.shadow.vtu                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)        :: grid
  character(*)           :: append
  real,         optional :: scalar_cell(-grid % n_bnd_cells:grid % n_cells)
  real,         optional :: scalar_node(1:grid % n_nodes)
  character(*), optional :: scalar_name
  real,         optional :: vector_cell(-grid % n_bnd_cells:grid % n_cells, 3)
  real,         optional :: vector_node(1:grid % n_nodes, 3)
  character(*), optional :: vector_name
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c, n, data_offset, cell_offset, n_conns, fu, lev
  character(SL) :: name_out, str1, str2
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: IP = DP  ! int. precision is double precision
  integer,           parameter :: RP = DP  ! real precision is double precision
  integer,           parameter :: VTK_LINE       =  3
  integer,           parameter :: VTK_TRIANGLE   =  5
  integer,           parameter :: VTK_QUAD       =  9
  integer,           parameter :: VTK_TETRA      = 10
  integer,           parameter :: VTK_HEXAHEDRON = 12
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
  character(len= 1), parameter :: LF   = char(10)      ! line feed
  character(len= 0), parameter :: IN_0 = ''            ! indentation levels
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '
!==============================================================================!

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do c = 1, grid % n_cells
    n_conns = n_conns + grid % cells_n_nodes(c)
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File_Mod_Set_Name(name_out, appendix='-'//trim(append),  &
                         processor=this_proc, extension='.vtu')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>'             // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid"'  //  &
                    ' version="0.1"'                    //  &
                    ' byte_order="LittleEndian">'       // LF
  write(fu) IN_1 // '<UnstructuredGrid>' // LF
  write(str1, '(i0.0)') grid % n_nodes
  write(str2, '(i0.0)') grid % n_cells
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"' //  &
                    ' NumberOfCells="' // trim(str2) // '">'       // LF
  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + grid % n_nodes * RP * 3  ! prepare for next

  !-----------!
  !   Cells   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Cells' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_cells * IP  ! prepare for next

  ! Cells' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_cells * IP  ! prepare for next

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !----------------!
  !   Point data   !
  !----------------!
  if(present(scalar_node) .or. present(vector_node)) then
    write(fu) IN_3 // '<PointData Scalars="scalars" vectors="velocity">' // LF
  end if

  ! Additional node-based scalar array
  if(present(scalar_node)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Float64"'            //  &
                      ' Name="'// trim(scalar_name) // '"'   //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + grid % n_nodes * RP  ! prepare for next
  end if

  ! Additional node-based vector array
  if(present(vector_node)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Float64"'            //  &
                      ' Name="'// trim(vector_name) // '"'   //  &
                      ' NumberOfComponents="3"'              //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + grid % n_nodes * RP * 3  ! prepare for next
  end if

  if(present(scalar_node) .or. present(vector_node)) then
    write(fu) IN_3 // '</PointData>' // LF
  end if

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  ! Processor i.d.
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="Processor"'              //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_cells * IP  ! prepare for next

  ! Additional cell scalar
  if(present(scalar_cell)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Float64"'            //  &
                      ' Name="'// trim(scalar_name) // '"'   //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + grid % n_cells * RP  ! prepare for next
  end if

  ! Additional cell vector
  if(present(vector_cell)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Float64"'            //  &
                      ' Name="'// trim(vector_name) // '"'   //  &
                      ' NumberOfComponents="3"'              //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + grid % n_cells * RP * 3  ! prepare for next
  end if

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu) IN_3 // '</CellData>'         // LF
  write(fu) IN_2 // '</Piece>'            // LF
  write(fu) IN_1 // '</UnstructuredGrid>' // LF

  !-------------------!
  !                   !
  !   Appended data   !
  !                   !
  !-------------------!
  write(fu) IN_0 // '<AppendedData encoding="raw">' // LF
  write(fu) '_'

  !-----------!
  !   Nodes   !
  !-----------!
  data_size = grid % n_nodes * RP * 3
  write(fu) data_size
  do n = 1, grid % n_nodes
    write(fu) grid % xn(n), grid % yn(n), grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Cells' nodes
  data_size = n_conns * IP
  write(fu) data_size
  do c = 1, grid % n_cells

    ! Hexahedral
    if(grid % cells_n_nodes(c) .eq. 8) then
      write(fu)                 &
        grid % cells_n(1,c)-1,  &
        grid % cells_n(2,c)-1,  &
        grid % cells_n(4,c)-1,  &
        grid % cells_n(3,c)-1,  &
        grid % cells_n(5,c)-1,  &
        grid % cells_n(6,c)-1,  &
        grid % cells_n(8,c)-1,  &
        grid % cells_n(7,c)-1

    ! Wedge
    else if(grid % cells_n_nodes(c) .eq. 6) then
      write(fu)                 &
        grid % cells_n(1,c)-1,  &
        grid % cells_n(2,c)-1,  &
        grid % cells_n(3,c)-1,  &
        grid % cells_n(4,c)-1,  &
        grid % cells_n(5,c)-1,  &
        grid % cells_n(6,c)-1

    ! Tetrahedra
    else if(grid % cells_n_nodes(c) .eq. 4) then
      write(fu)                 &
        grid % cells_n(1,c)-1,  &
        grid % cells_n(2,c)-1,  &
        grid % cells_n(3,c)-1,  &
        grid % cells_n(4,c)-1

    ! Pyramid
    else if(grid % cells_n_nodes(c) .eq. 5) then
      write(fu)                 &
        grid % cells_n(1,c)-1,  &
        grid % cells_n(2,c)-1,  &
        grid % cells_n(4,c)-1,  &
        grid % cells_n(3,c)-1,  &
        grid % cells_n(5,c)-1
    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if
  end do

  ! Cells' offsets
  data_size = grid % n_cells * IP
  write(fu) data_size
  cell_offset = 0
  do c = 1, grid % n_cells
    cell_offset = cell_offset + grid % cells_n_nodes(c)
    write(fu) cell_offset
  end do

  ! Cells' types
  data_size = grid % n_cells * IP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) .eq. 4) write(fu) VTK_TETRA
    if(grid % cells_n_nodes(c) .eq. 8) write(fu) VTK_HEXAHEDRON
    if(grid % cells_n_nodes(c) .eq. 6) write(fu) VTK_WEDGE
    if(grid % cells_n_nodes(c) .eq. 5) write(fu) VTK_PYRAMID
  end do

  !----------------!
  !   Point data   !
  !----------------!
  if(present(scalar_node)) then
    data_size = grid % n_nodes * RP
    write(fu) data_size
    do n = 1, grid % n_nodes
      write(fu) scalar_node(n)
    end do
  end if

  if(present(vector_node)) then
    data_size = grid % n_nodes * RP * 3
    write(fu) data_size
    do n = 1, grid % n_nodes
      write(fu) vector_node(n, 1), vector_node(n, 2), vector_node(n, 3)
    end do
  end if

  !---------------!
  !   Cell data   !
  !---------------!

  ! Processor i.d.
  data_size = grid % n_cells * IP
  write(fu) data_size
  do c = 1, grid % n_cells
    write(fu) grid % comm % cell_proc(c)
  end do

  ! Additional cell data
  if(present(scalar_cell)) then
    data_size = grid % n_cells * RP
    write(fu) data_size
    do c = 1, grid % n_cells
      write(fu) scalar_cell(c)
    end do
  end if

  if(present(vector_cell)) then
    data_size = grid % n_cells * RP * 3
    write(fu) data_size
    do c = 1, grid % n_cells
      write(fu) vector_cell(c, 1), vector_cell(c, 2), vector_cell(c, 3)
    end do
  end if

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!

  ! Create it only from subdomain 1, when decomposed
  if(maxval(grid % comm % cell_proc(:)) > 1 .and. this_proc .eq. 1) then

    call File_Mod_Set_Name(name_out, appendix='-'//trim(append),  &
                           extension='.pvtu')
    call File_Mod_Open_File_For_Writing(name_out, fu)

    ! Header
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(fu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(fu,'(a,a)') IN_2, '<PPoints>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" NumberOfComponents='// &
                           '"3" format="ascii"/>'
    write(fu,'(a,a)') IN_2, '</PPoints>'

    ! Data section is not mandatory, but very useful
    write(fu,'(a,a)') IN_2, '<PCellData Scalars="scalars" vectors="velocity">'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Int64" Name="Processor"' // &
                           ' format="ascii"/>'
    if(present(scalar_cell)) then
      write(fu,'(a,a)') IN_3, '<PDataArray type="Float64"'         //  &
                              ' Name="'// trim(scalar_name) // '"' //  &
                              ' format="ascii"/>'
    end if
    if(present(vector_cell)) then
      write(fu,'(a,a)') IN_3, '<PDataArray type="Float64"'         //  &
                              ' NumberOfComponents="3"'            //  &
                              ' Name="'// trim(vector_name) // '"' //  &
                              ' format="ascii"/>'
    end if
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, n_proc
      call File_Mod_Set_Name(name_out, appendix='-'//trim(append),  &
                             processor=n, extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  end subroutine
