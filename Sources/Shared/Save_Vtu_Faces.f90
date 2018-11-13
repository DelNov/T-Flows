!==============================================================================!
  subroutine Save_Vtu_Faces(grid)
!------------------------------------------------------------------------------!
! Writes .faces.vtu file.                                                      !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c1, c2, n, s, offset
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
  call Name_File(0, name_out, '.faces.vtu')
  open(9, file=name_out)
  print *, '# Creating the file: ', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                         'byte_order="LittleEndian">'
  write(9,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(9,'(a,a,i0.0,a,i0.0,a)')   &
                   IN_2, '<Piece NumberOfPoints="', grid % n_nodes, &
                              '" NumberOfCells ="', grid % n_faces, '">'
  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Points>'
  write(9,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents=' //  &
                 '"3" format="ascii">'
  do n = 1, grid % n_nodes
    write(9, '(a,1PE15.7,1PE15.7,1PE15.7)')                &
               IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Faces   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Cells>'

  ! First write all faces' nodes
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                         ' format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) .eq. 4) then
      write(9,'(a,4I9)')                               &
        IN_5,                                          &
        grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
        grid % faces_n(3,s)-1, grid % faces_n(4,s)-1
    else if(grid % faces_n_nodes(s) .eq. 3) then
      write(9,'(a,3I9)')                               &
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
  write(9,'(a,a)') IN_4, '</DataArray>'

  ! Then write all faces' offsets
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets" format="ascii">'
  offset = 0
  do s = 1, grid % n_faces
    offset = offset + grid % faces_n_nodes(s)
    write(9,'(a,i9)') IN_5, offset
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) .eq. 4) write(9,'(a,i9)') IN_5, VTK_QUAD
    if(grid % faces_n_nodes(s) .eq. 3) write(9,'(a,i9)') IN_5, VTK_TRIANGLE
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !   Cell data   !
  !---------------!
  write(9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  ! Boundary conditions
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" ' // &
                   'Name="BoundaryConditions" format="ascii">'
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! If boundary
    if( c2 < 0 ) then
      write(9,'(a,i9)') IN_5, grid % bnd_cond % color(c2)

    ! If inside
    else
      write(9,'(a,i9)') IN_5, 0
    end if
  end do

  write(9,'(a,a)') IN_4, '</DataArray>'

  !------------!
  !   Footer   !
  !------------!
  write(9,'(a,a)') IN_3, '</CellData>'
  write(9,'(a,a)') IN_2, '</Piece>'
  write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(9,'(a,a)') IN_0, '</VTKFile>'

  close(9)

  end subroutine
