!==============================================================================!
  subroutine Save_Poly_Vtu(grid)
!------------------------------------------------------------------------------!
!   Writes: name.vtu, name.faces.vtu, name.shadow.vtu                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, i_pol, offset, fu
  character(SL) :: name_out, str1, str2, str3
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: VTK_TETRA      = 10  ! cells in VTK format
  integer,           parameter :: VTK_HEXAHEDRON = 12
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
  integer,           parameter :: VTK_POLYHEDRON = 42
  character(len= 1), parameter :: LF   = char(10)      ! line feed
  character(len= 0), parameter :: IN_0 = ''            ! indentation levels 
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, extension='.vtu')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !------------!
  !   Header   !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>'             // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid" ' //  &
                    'version="0.1" '                    //  &
                    'byte_order="LittleEndian">'        // LF
  write(fu) IN_1 // '<UnstructuredGrid>'                // LF
  write(str1, '(i0.0)') grid % n_nodes
  write(str2, '(i0.0)') grid % n_cells
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"' //  &
                    ' NumberOfCells="' // trim(str2) // '">'       // LF

  !-----------!
  !   Nodes   !
  !-----------!
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="ascii">'               // LF
  do n = 1, grid % n_nodes
    write(str1, '(1pe15.7)') grid % xn(n)
    write(str2, '(1pe15.7)') grid % yn(n)
    write(str3, '(1pe15.7)') grid % zn(n)
    write(fu) IN_5 // trim(str1) // trim(str2) // trim(str3) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF

  !-----------!
  !   Cells   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! First write all cells' nodes
  write(fu) IN_4 // '<DataArray type="Int64"' //  &
                    ' Name="connectivity"'    //  &
                    ' format="ascii">'        // LF

  do c = 1, grid % n_cells

    ! Hexahedral
    if(grid % cells_n_nodes(c) .eq. 8) then
      write(str1,'(64i9)') (grid % cells_n(1:grid % cells_n_nodes(c), c))-1
      write(fu) IN_5 // trim(str1) // LF

    ! Wedge
    else if(grid % cells_n_nodes(c) .eq. 6) then
      write(str1,'(64i9)') (grid % cells_n(1:grid % cells_n_nodes(c), c))-1
      write(fu) IN_5 // trim(str1) // LF

    ! Tetrahedra
    else if(grid % cells_n_nodes(c) .eq. 4) then
      write(str1,'(64i9)') (grid % cells_n(1:grid % cells_n_nodes(c), c))-1
      write(fu) IN_5 // trim(str1) // LF

    ! Pyramid
    else if(grid % cells_n_nodes(c) .eq. 5) then
      write(str1,'(64i9)') (grid % cells_n(1:grid % cells_n_nodes(c), c))-1
      write(fu) IN_5 // trim(str1) // LF

    ! Polyhedral cells
    else if(grid % cells_n_nodes(c) < 0) then
      write(str1,'(64i9)') (grid % cells_n(1:grid % cells_n_nodes(c), c))-1
      write(fu) IN_5 // trim(str1) // LF

    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop
    end if

  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Now write all cells' offsets
  write(fu) IN_4 // '<DataArray type="Int64"'  //  &
                    ' Name="offsets"'          //  &
                    ' format="ascii">'         // LF
  offset = 0
  do c = 1, grid % n_cells
    offset = offset + abs(grid % cells_n_nodes(c))
    write(str1,'(i9)') offset
    write(fu) IN_5 // trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Now write all cells' types
  write(fu) IN_4 // '<DataArray type="Int64"'  //  &
                    ' Name="types"'            //  &
                    ' format="ascii">'         // LF
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) .eq. 4) write(str1,'(i9)') VTK_TETRA
    if(grid % cells_n_nodes(c) .eq. 8) write(str1,'(i9)') VTK_HEXAHEDRON
    if(grid % cells_n_nodes(c) .eq. 6) write(str1,'(i9)') VTK_WEDGE
    if(grid % cells_n_nodes(c) .eq. 5) write(str1,'(i9)') VTK_PYRAMID
    if(grid % cells_n_nodes(c) .lt. 0) write(str1,'(i9)') VTK_POLYHEDRON
    write(fu) IN_5 // trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Write polyhedral cells' faces
  write(fu) IN_4 // '<DataArray type="Int64"'  //  &
                    ' Name="faces"'            //  &
                    ' format="ascii">'         // LF
  do c = 1, grid % n_cells

    ! You have found a polyhedron, write its faces out
    if(grid % cells_n_nodes(c) .lt. 0) then

      ! Write number of polyfaces for this cell
      write(str1,'(i9)') grid % cells_n_polyf(c)
      write(fu) IN_5 // trim(str1) // LF

      do i_pol = 1, grid % cells_n_polyf(c)
        s = grid % cells_p(i_pol, c)
        n = grid % faces_n_nodes(s)
        write(str1,'(64i9)')  grid % faces_n_nodes(s),  &
                             (grid % faces_n(1:n, s))-1
        write(fu) IN_5 // trim(str1) // LF
      end do
    end if
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Write polyhedral cells' faces offsets
  offset = 0
  write(fu) IN_4 // '<DataArray type="Int64"'  //  &
                    ' Name="facesoffsets"'     //  &
                    ' format="ascii">'         // LF
  do c = 1, grid % n_cells

    ! You have found a polyhedron
    if(grid % cells_n_nodes(c) .lt. 0) then

      ! Increase offset for storing number of polyfaces
      offset = offset + 1

      ! Update the offset with all faces and their nodes
      do i_pol = 1, grid % cells_n_polyf(c)
        s = grid % cells_p(i_pol, c)
        n = grid % faces_n_nodes(s)
        offset = offset + 1 + n
      end do

      ! Write the current offset
      write(str1,'(i9)') offset
      write(fu) IN_5 // trim(str1) // LF

    ! Not a polyhedron, offsets are not needed
    else
      write(str1,'(i9)') -1
      write(fu) IN_5 // trim(str1) // LF
    end if

  end do
  write(fu) IN_4 // '</DataArray>' // LF

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  ! Processor i.d.
  write(fu) IN_4 // '<DataArray type="Int64"'           //  &
                    ' Name="Processor" format="ascii">' // LF
  do c = 1, grid % n_cells
    write(str1,'(i9)') grid % comm % cell_proc(c)
    write(fu) IN_5 // trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Wall distance
  write(fu) IN_4 // '<DataArray type="Float64"'            //  &
                    ' Name="WallDistance" format="ascii">' // LF
  do c = 1, grid % n_cells
    write(str1,'(1pe15.7)') grid % wall_dist(c)
    write(fu) IN_5 // trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  ! Cell volume
  write(fu) IN_4 // '<DataArray type="Float64"'          //  &
                    ' Name="CellVolume" format="ascii">' // LF
  do c = 1, grid % n_cells
    write(str1,'(1pe15.7)') grid % vol(c)
    write(fu) IN_5 // trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF

  !------------!
  !   Footer   !
  !------------!
  write(fu) IN_3 // '</CellData>'         // LF
  write(fu) IN_2 // '</Piece>'            // LF
  write(fu) IN_1 // '</UnstructuredGrid>' // LF
  write(fu) IN_0 // '</VTKFile>'          // LF

  close(fu)

  end subroutine
