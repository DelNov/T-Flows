!==============================================================================!
  subroutine Save_Vtu_Faces(grid)
!------------------------------------------------------------------------------!
!   Writes .faces.vtu file.                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c1, c2, n, s, offset, fu
  character(len=80)  :: name_out
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: VTK_TRIANGLE = 5  ! cell shapes in VTK format
  integer,           parameter :: VTK_QUAD     = 9
  character(len= 0), parameter :: IN_0 = ''         ! indentation levels 
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '
!==============================================================================!

  !-----------------------------------------!
  !                                         !
  !   Create boundary condition .vtu file   !
  !                                         !
  !-----------------------------------------!
  call File_Mod_Set_Name(name_out, extension='.faces.vtu')
  call File_Mod_Open_File_For_Writing(name_out, fu)

  !-----------!
  !   Start   !
  !-----------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                         'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                    IN_2, '<Piece NumberOfPoints="', grid % n_nodes, &
                               '" NumberOfCells ="', grid % n_faces, '">'
  !-----------!
  !   Nodes   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents=' //  &
                 '"3" format="ascii">'
  do n = 1, grid % n_nodes
    write(fu, '(a,1PE15.7,1PE15.7,1PE15.7)')                        &
                    IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Faces   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'

  ! First write all faces' nodes
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                         ' format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) .eq. 4) then
      write(fu,'(a,4I9)')                              &
        IN_5,                                          &
        grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
        grid % faces_n(3,s)-1, grid % faces_n(4,s)-1
    else if(grid % faces_n_nodes(s) .eq. 3) then
      write(fu,'(a,3I9)')                              &
        IN_5,                                          &
        grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
        grid % faces_n(3,s)-1
    else
      print *, '# Unsupported cell type ',       &
                 grid % faces_n_nodes(s), ' nodes.'
      print *, '# Exiting'
      stop
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Then write all faces' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets" ' // &
                          'format="ascii">'
  offset = 0
  do s = 1, grid % n_faces
    offset = offset + grid % faces_n_nodes(s)
    write(fu,'(a,i9)') IN_5, offset
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_QUAD
    if(grid % faces_n_nodes(s) .eq. 3) write(fu,'(a,i9)') IN_5, VTK_TRIANGLE
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  ! Boundary conditions
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' // &
                    'Name="BoundaryConditions" format="ascii">'
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! If boundary
    if( c2 < 0 ) then
      write(fu,'(a,i9)') IN_5, grid % bnd_cond % color(c2)

    ! If inside
    else
      write(fu,'(a,i9)') IN_5, 0
    end if
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
