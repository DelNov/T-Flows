!==============================================================================!
  subroutine Save_Vtu_Edges(grid, edge_data)
!------------------------------------------------------------------------------!
!   Writes file with edges                                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  integer, optional :: edge_data(grid % n_edges)
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, edge_offset, fu
  character(SL) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, extension='.edges.vtu')
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
                    IN_2, '<Piece NumberOfPoints="', grid % n_nodes,      &
                               '" NumberOfCells ="', grid % n_edges, '">'

  !-----------!
  !           !
  !   Nodes   !
  !           !
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
  !           !
  !   Edges   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'

  ! First write all edges' nodes
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                          ' format="ascii">'
  do c = 1, grid % n_edges
    write(fu,'(a,64i9)') IN_5, (grid % edges_n(1:2, c))-1
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all edges' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                          'Name="offsets" format="ascii">'
  edge_offset = 0
  do c = 1, grid % n_edges
    edge_offset = edge_offset + 2
    write(fu,'(a,i9)') IN_5, edge_offset
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all edges' types
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do c = 1, grid % n_edges
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !               !
  !   Edge data   !
  !               !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  if(present(edge_data)) then
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                            'Name="Processor" format="ascii">'
    do c = 1, grid % n_edges
      write(fu,'(a,i9)') IN_5, edge_data(c)
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
  end if
  !Not yet:
  !Not yet:  ! Wall distance
  !Not yet:  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
  !Not yet:                          'Name="GeomWallDistance" format="ascii">'
  !Not yet:  do c = 1, grid % n_cells
  !Not yet:    write(fu,'(a,1pe15.7)') IN_5, grid % wall_dist(c)
  !Not yet:  end do
  !Not yet:  write(fu,'(a,a)') IN_4, '</DataArray>'
  !Not yet:
  !Not yet:  ! Cell volume
  !Not yet:  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
  !Not yet:                          'Name="GeomCellVolume" format="ascii">'
  !Not yet:  do c = 1, grid % n_cells
  !Not yet:    write(fu,'(a,1pe15.7)') IN_5, grid % vol(c)
  !Not yet:  end do
  !Not yet:  write(fu,'(a,a)') IN_4, '</DataArray>'

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
