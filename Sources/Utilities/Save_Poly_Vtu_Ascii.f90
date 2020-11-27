!==============================================================================!
  subroutine Save_Poly_Vtu_Ascii(grid)
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
  character(SL) :: name_out
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: VTK_TETRA      = 10  ! cells in VTK format
  integer,           parameter :: VTK_HEXAHEDRON = 12
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
  integer,           parameter :: VTK_POLYHEDRON = 42
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
  call File_Mod_Open_File_For_Writing(name_out, fu)

  !------------!
  !   Header   !
  !------------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//  &
                          'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                    IN_2, '<Piece NumberOfPoints="', grid % n_nodes,      &
                               '" NumberOfCells ="', grid % n_cells, '">'

  !-----------!
  !   Nodes   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                          '="3" format="ascii">'
  do n = 1, grid % n_nodes
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')  &
              IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Cells   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'

  ! First write all cells' nodes
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                          ' format="ascii">'

  do c = 1, grid % n_cells

    ! Hexahedral
    if(grid % cells_n_nodes(c) .eq. 8) then
      write(fu,'(a,64i9)')  &
        IN_5, (grid % cells_n(1:grid % cells_n_nodes(c), c))-1

    ! Wedge
    else if(grid % cells_n_nodes(c) .eq. 6) then
      write(fu,'(a,64i9)')  &
        IN_5, (grid % cells_n(1:grid % cells_n_nodes(c), c))-1

    ! Tetrahedra
    else if(grid % cells_n_nodes(c) .eq. 4) then
      write(fu,'(a,64i9)')  &
        IN_5, (grid % cells_n(1:grid % cells_n_nodes(c), c))-1

    ! Pyramid
    else if(grid % cells_n_nodes(c) .eq. 5) then
      write(fu,'(a,64i9)')  &
        IN_5, (grid % cells_n(1:grid % cells_n_nodes(c), c))-1

    ! Polyhedral cells
    else if(grid % cells_n_nodes(c) < 0) then
      write(fu,'(a,64i9)')  &
        IN_5, (grid % cells_n(1:-grid % cells_n_nodes(c), c))-1

    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if

  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                          'Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    offset = offset + abs(grid % cells_n_nodes(c))
    write(fu,'(a,i9)') IN_5, offset
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_TETRA
    if(grid % cells_n_nodes(c) .eq. 8) write(fu,'(a,i9)') IN_5, VTK_HEXAHEDRON
    if(grid % cells_n_nodes(c) .eq. 6) write(fu,'(a,i9)') IN_5, VTK_WEDGE
    if(grid % cells_n_nodes(c) .eq. 5) write(fu,'(a,i9)') IN_5, VTK_PYRAMID
    if(grid % cells_n_nodes(c) .lt. 0) write(fu,'(a,i9)') IN_5, VTK_POLYHEDRON
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Write polyhedral cells' faces
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="faces" format="ascii">'
  do c = 1, grid % n_cells

    ! You have found a polyhedron, write its faces out
    if(grid % cells_n_nodes(c) .lt. 0) then

      ! Write number of polyfaces for this cell
      write(fu,'(a,i9)') IN_5, grid % cells_n_polyf(c)

      do i_pol = 1, grid % cells_n_polyf(c)
        s = grid % cells_p(i_pol, c)
        n = grid % faces_n_nodes(s)
        write(fu,'(a,64i9)') IN_5,  grid % faces_n_nodes(s),  &
                                   (grid % faces_n(1:n, s))-1
      end do
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Write polyhedral cells' faces offsets
  offset = 0
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="faceoffsets" format="ascii">'
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
      write(fu,'(a,i9)') IN_5, offset

    ! Not a polyhedron, offsets are not needed
    else
      write(fu,'(a,i9)') IN_5, -1
    end if

  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  ! Processor i.d.
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                          'Name="Processor" format="ascii">'
  do c = 1, grid % n_cells
    write(fu,'(a,i9)') IN_5, grid % comm % cell_proc(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Wall distance
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
                          'Name="WallDistance" format="ascii">'
  do c = 1, grid % n_cells
    write(fu,'(a,1pe15.7)') IN_5, grid % wall_dist(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Cell volume
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
                          'Name="CellVolume" format="ascii">'
  do c = 1, grid % n_cells
    write(fu,'(a,1pe15.7)') IN_5, grid % vol(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !------------!
  !   Footer   !
  !------------!
  write(fu,'(a,a)') IN_3, '</CellData>'
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'

  close(fu)

  end subroutine
