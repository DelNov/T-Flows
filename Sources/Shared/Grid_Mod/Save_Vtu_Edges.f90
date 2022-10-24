!==============================================================================!
  subroutine Save_Vtu_Edges(Grid, edge_data)
!------------------------------------------------------------------------------!
!   Writes edges in vtu file format                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)  :: Grid
  integer, optional :: edge_data(Grid % n_edges)
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, edge_offset, fu
  character(SL) :: name_out
!==============================================================================!

  call Profiler % Start('Save_Vtu_Edges')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call File % Set_Name(name_out, extension='.edges.vtu')
  call File % Open_For_Writing_Ascii(name_out, fu)

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
                               '" NumberOfCells ="', Grid % n_edges, '">'

  !-----------!
  !           !
  !   Nodes   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp//' NumberOfComponents' // &
                          '="3" format="ascii">'
  do n = 1, Grid % n_nodes
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')  &
              IN_5, Grid % xn(n), Grid % yn(n), Grid % zn(n)
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
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//  &
                          ' Name="connectivity"' //  &
                          ' format="ascii">'
  do c = 1, Grid % n_edges
    write(fu,'(a,64i9)') IN_5, (Grid % edges_n(1:2, c))-1
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all edges' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//  &
                          ' Name="offsets" format="ascii">'
  edge_offset = 0
  do c = 1, Grid % n_edges
    edge_offset = edge_offset + 2
    write(fu,'(a,i9)') IN_5, edge_offset
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all edges' types
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//  &
                         ' Name="types" format="ascii">'
  do c = 1, Grid % n_edges
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !               !
  !   Edge data   !
  !               !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData>'

  if(present(edge_data)) then
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//  &
                            ' Name="EdgeData" format="ascii">'
    do c = 1, Grid % n_edges
      write(fu,'(a,i9)') IN_5, edge_data(c)
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
  end if

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

  call Profiler % Stop('Save_Vtu_Edges')

  end subroutine
