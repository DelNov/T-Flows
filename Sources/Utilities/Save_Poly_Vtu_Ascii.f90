!==============================================================================!
  subroutine Save_Poly_Vtu_Ascii(Grid)
!------------------------------------------------------------------------------!
!   Writes: name.vtu, name.faces.vtu, name.shadow.vtu                          !
!                                                                              !
!   If you change precision of integers to 64 bits (currently set in the       !
!   makefile and in Const_Mod.f90 with parameter IP), all occurences of        !
!   Int32 here should be changed to Int64.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, i_pol, cell_offset, fu
  character(SL) :: name_out
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, extension='.vtu')
  call File_Mod_Open_File_For_Writing(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//  &
                          'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                    IN_2, '<Piece NumberOfPoints="', Grid % n_nodes,      &
                               '" NumberOfCells ="', Grid % n_cells, '">'

  !-----------!
  !           !
  !   Nodes   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          ' NumberOfComponents=3'     //  &
                          ' format="ascii">'
  do n = 1, Grid % n_nodes
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')  &
              IN_5, Grid % xn(n), Grid % yn(n), Grid % zn(n)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'

  ! First write all cells' nodes
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="connectivity"'    //  &
                          ' format="ascii">'

  do c = 1, Grid % n_cells

    ! Hexahedral
    if(Grid % cells_n_nodes(c) .eq. 8) then
      write(fu,'(a,64i9)')  &
        IN_5, (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

    ! Wedge
    else if(Grid % cells_n_nodes(c) .eq. 6) then
      write(fu,'(a,64i9)')  &
        IN_5, (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

    ! Tetrahedra
    else if(Grid % cells_n_nodes(c) .eq. 4) then
      write(fu,'(a,64i9)')  &
        IN_5, (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

    ! Pyramid
    else if(Grid % cells_n_nodes(c) .eq. 5) then
      write(fu,'(a,64i9)')  &
        IN_5, (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

    ! Polyhedral cells
    else if(Grid % cells_n_nodes(c) < 0) then
      write(fu,'(a,64i9)')  &
        IN_5, (Grid % cells_n(1:-Grid % cells_n_nodes(c), c))-1

    else
      print *, '# Unsupported cell type with ',  &
                  Grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if

  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          'Name="offsets" format="ascii">'
  cell_offset = 0
  do c = 1, Grid % n_cells
    cell_offset = cell_offset + abs(Grid % cells_n_nodes(c))
    write(fu,'(a,i9)') IN_5, cell_offset
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="types" format="ascii">'
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_TETRA
    if(Grid % cells_n_nodes(c) .eq. 8) write(fu,'(a,i9)') IN_5, VTK_HEXAHEDRON
    if(Grid % cells_n_nodes(c) .eq. 6) write(fu,'(a,i9)') IN_5, VTK_WEDGE
    if(Grid % cells_n_nodes(c) .eq. 5) write(fu,'(a,i9)') IN_5, VTK_PYRAMID
    if(Grid % cells_n_nodes(c) .lt. 0) write(fu,'(a,i9)') IN_5, VTK_POLYHEDRON
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Write polyhedral cells' faces
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="faces" format="ascii">'
  do c = 1, Grid % n_cells

    ! You have found a polyhedron, write its faces out
    if(Grid % cells_n_nodes(c) .lt. 0) then

      ! Write number of polyfaces for this cell
      write(fu,'(a,i9)') IN_5, Grid % cells_n_faces(c)

      do i_pol = 1, Grid % cells_n_faces(c)
        s = Grid % cells_f(i_pol, c)
        n = Grid % faces_n_nodes(s)
        write(fu,'(a,64i9)') IN_5,  Grid % faces_n_nodes(s),  &
                                   (Grid % faces_n(1:n, s))-1
      end do
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Write polyhedral cells' faces offsets
  cell_offset = 0
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="faceoffsets" format="ascii">'
  do c = 1, Grid % n_cells

    ! You have found a polyhedron
    if(Grid % cells_n_nodes(c) .lt. 0) then

      ! Increase offset for storing number of polyfaces
      cell_offset = cell_offset + 1

      ! Update the offset with all faces and their nodes
      do i_pol = 1, Grid % cells_n_faces(c)
        s = Grid % cells_f(i_pol, c)
        n = Grid % faces_n_nodes(s)
        cell_offset = cell_offset + 1 + n
      end do

      ! Write the current offset
      write(fu,'(a,i9)') IN_5, cell_offset

    ! Not a polyhedron, offsets are not needed
    else
      write(fu,'(a,i9)') IN_5, -1
    end if

  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData>'

  ! Processor i.d.
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          'Name="Processor" format="ascii">'
  do c = 1, Grid % n_cells
    write(fu,'(a,i9)') IN_5, Grid % Comm % cell_proc(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Wall distance
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          'Name="GeomWallDistance" format="ascii">'
  do c = 1, Grid % n_cells
    write(fu,'(a,1pe15.7)') IN_5, Grid % wall_dist(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Cell volume
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          'Name="GeomCellVolume" format="ascii">'
  do c = 1, Grid % n_cells
    write(fu,'(a,1pe15.7)') IN_5, Grid % vol(c)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  write(fu,'(a,a)') IN_3, '</CellData>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'

  !---------------------!
  !                     !
  !   Close .vtu file   !
  !                     !
  !---------------------!
  close(fu)

  end subroutine
