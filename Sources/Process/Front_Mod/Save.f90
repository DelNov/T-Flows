!==============================================================================!
  subroutine Front_Mod_Save(front, time_step)
!------------------------------------------------------------------------------!
!   Writes surface vertices in VTU file format (for VisIt and Paraview)        !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Front_Type), target :: front
  integer                  :: time_step
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: vert
  integer                  :: v, e     ! vertex and element counters
  integer                  :: offset, fu
  character(SL)            :: name_out
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: VTK_TRIANGLE = 5  ! cell shapes in VTK format
  character(len= 0)  :: IN_0 = ''         ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(front % n_verts < 1) return

  !----------------------------!
  !                            !
  !   Create .front.vtu file   !
  !                            !
  !----------------------------!

  if(this_proc < 2) then

    call File_Mod_Set_Name(name_out,               &
                           time_step = time_step,  &
                           appendix  = '-front',   &
                           extension = '.vtu')
    call File_Mod_Open_File_For_Writing(name_out, fu)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                            'byte_order="LittleEndian">'
    write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                IN_2, '<Piece NumberOfPoints="', front % n_verts,  &
                           '" NumberOfCells ="', front % n_elems, '">'

    !------------------------!
    !                        !
    !   Vertex coordinates   !
    !                        !
    !------------------------!
    write(fu,'(a,a)') IN_3, '<Points>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                            '="3" format="ascii">'
    do v = 1, front % n_verts
      vert => front % vert(v)
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                  IN_5, vert % x_n, vert % y_n, vert % z_n
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Points>'

    !----------------!
    !                !
    !   Point data   !
    !                !
    !----------------!
    write(fu,'(a,a)') IN_3, '<PointData Scalars="scalars" vectors="velocity">'

    !--------------------!
    !   Particle i.d.s   !
    !--------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Index" ' // &
                            'format="ascii">'
    do v = 1, front % n_verts
      write(fu,'(a,i9)') IN_5, v
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !--------------------------!
    !   Number of neighbours   !
    !--------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours" ' // &
                            'format="ascii">'
    do v = 1, front % n_verts
      write(fu,'(a,i9)') IN_5, front % vert(v) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-----------------------------!
    !   Curvatures at the nodes   !
    !-----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" Name="NodeCurv" ' // &
                           ' format="ascii">'
    do v = 1, front % n_verts
      vert => front % vert(v)
      write(fu,'(a,1pe16.6e4)') IN_5, vert % curv
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    write(fu,'(a,a)') IN_3, '</PointData>'

    !-----------!
    !           !
    !   Cells   !
    !           !
    !-----------!
    write(fu,'(a,a)') IN_3, '<Cells>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                            ' format="ascii">'
    ! Cell topology
    do e = 1, front % n_elems
      write(fu,'(a,99i9)') IN_5, front % elem(e) % v(1:front % elem(e) % nv)-1
    end do

    ! Cell offsets
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                            ' format="ascii">'
    offset = 0
    do e = 1, front % n_elems
      offset = offset + front % elem(e) % nv
      write(fu,'(a,i9)') IN_5, offset
    end do

    ! Cell types
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types"' //  &
                            ' format="ascii">'
    do e = 1, front % n_elems
      write(fu,'(a,i9)') IN_5, VTK_POLYGON
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Cells>'

    !---------------!
    !               !
    !   Cell data   !
    !               !
    !---------------!

    ! Beginning of cell data
    write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

    !-------------------------------------!
    !   Number of neighbouring elements   !
    !-------------------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours"' //  &
                            ' format="ascii">'
    do e = 1, front % n_elems
      write(fu,'(a,i9)') IN_5, front % elem(e) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !---------------------!
    !   Surface normals   !
    !---------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementNormals" ' //  &
                   ' NumberOfComponents="3" format="ascii">'
    do e = 1, front % n_elems
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                IN_5, front % elem(e) % nx,           &
                      front % elem(e) % ny,           &
                      front % elem(e) % nz
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------!
    !   Element areas   !
    !-------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementArea" ' //  &
                   ' format="ascii">'
    do e = 1, front % n_elems
      write(fu,'(a,1pe16.6e4)') IN_5, front % elem(e) % area
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Element coordinates   !
    !-------------------------!
    write(fu,'(4a)') IN_4,                                                    &
                   '<DataArray type="Float64" Name="ElementCoordinates" ' //  &
                   ' NumberOfComponents="3" format="ascii">'
    do e = 1, front % n_elems
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                IN_5, front % elem(e) % xe,           &
                      front % elem(e) % ye,           &
                      front % elem(e) % ze
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------!
    !   Surface curvatures   !
    !------------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementCurv" ' //  &
                   ' format="ascii">'
    do e = 1, front % n_elems
      write(fu,'(a,1pe16.6e4)') IN_5, front % elem(e) % curv
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    ! End of cell data
    write(fu,'(a,a)') IN_3, '</CellData>'

    !------------!
    !            !
    !   Footer   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_2, '</Piece>'
    write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
    write(fu,'(a,a)') IN_0, '</VTKFile>'
    close(fu)
  end if

  end subroutine
