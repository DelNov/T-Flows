!==============================================================================!
  subroutine Save_Surf(surf, time_step)
!------------------------------------------------------------------------------!
!   Writes surface vertices in VTU file format (for VisIt and Paraview)        !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Surf_Type), target :: surf
  integer                 :: time_step
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: vert
  integer                  :: v, e     ! vertex and element counters
  integer                  :: offset
  character(len=80)        :: name_out_9
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: VTK_TRIANGLE = 5  ! cell shapes in VTK format
  character(len= 0)  :: IN_0 = ''         ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(surf % n_verts < 1) return

  !---------------------------!
  !                           !
  !   Create .surf.vtu file   !
  !                           !
  !---------------------------!

  if(this_proc < 2) then

    call File_Mod_Set_Name(name_out_9, time_step=time_step, extension='.surf.vtu')

    open(9, file=name_out_9)
    print *, '# Creating file: ', trim(name_out_9)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //&
                           'byte_order="LittleEndian">'
    write(9,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(9,'(a,a,i0.0,a,i0.0,a)')   &
               IN_2, '<Piece NumberOfPoints="', surf % n_verts,  &
                          '" NumberOfCells ="', surf % n_elems, '">'

    !------------------------!
    !                        !
    !   Vertex coordinates   !
    !                        !
    !------------------------!
    write(9,'(a,a)') IN_3, '<Points>'
    write(9,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                           '="3" format="ascii">'
    do v = 1, surf % n_verts
      vert => surf % vert(v)
      write(9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                 IN_5, vert % x_n, vert % y_n, vert % z_n
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_3, '</Points>'

    !----------------!
    !                !
    !   Point data   !
    !                !
    !----------------!
    write(9,'(a,a)') IN_3, '<PointData Scalars="scalars" vectors="velocity">'

    !--------------------!
    !   Particle i.d.s   !
    !--------------------!
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="Index" ' // &
                           'format="ascii">'
    do v = 1, surf % n_verts
      write(9,'(a,i9)') IN_5, v
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'

    !--------------------------!
    !   Number of neighbours   !
    !--------------------------!
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours" ' // &
                           'format="ascii">'
    do v = 1, surf % n_verts
      write(9,'(a,i9)') IN_5, surf % vert(v) % nne
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'

!   !-------------------------!
!   !   Particle velocities   !
!   !-------------------------!
!   write(9,'(a,a)') IN_4, '<DataArray type="Float64" Name="Velocity" ' // &
!                          ' NumberOfComponents="3" format="ascii">'
!   do v = 1, surf % n_verts
!     part => surf % vert(k)
!     write(9,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                         &
!               IN_5, part % u, part % v, part % w
!   end do
!   write(9,'(a,a)') IN_4, '</DataArray>'
!
    write(9,'(a,a)') IN_3, '</PointData>'

    !-----------!
    !           !
    !   Cells   !
    !           !
    !-----------!
    write(9,'(a,a)') IN_3, '<Cells>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                           ' format="ascii">'
    ! Cell topology
    do e = 1, surf % n_elems
      write(9,'(a,3i9)')       &
        IN_5,                  &
        surf % elem(e) % i-1,  &
        surf % elem(e) % j-1,  &
        surf % elem(e) % k-1
    end do

    ! Cell offsets
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                           ' format="ascii">'
    offset = 0
    do e = 1, surf % n_elems
      offset = offset + 3
      write(9,'(a,i9)') IN_5, offset
    end do

    ! Cell types
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="types"' //  &
                           ' format="ascii">'
    do e = 1, surf % n_elems
      write(9,'(a,i9)') IN_5, VTK_TRIANGLE
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_3, '</Cells>'

    !---------------!
    !               !
    !   Cell data   !
    !               !
    !---------------!

    ! Beginning of cell data
    write(9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

    !-------------------------------------!
    !   Number of neighbouring elements   !
    !-------------------------------------!
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours"' //  &
                           ' format="ascii">'
    do e = 1, surf % n_elems
      write(9,'(a,i9)') IN_5, surf % elem(e) % nne
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'

    !---------------------!
    !   Surface normals   !
    !---------------------!
    write(9,'(4a)') IN_4,                                                &
                  '<DataArray type="Float64" Name="SurfaceNormals" ' //  &
                  ' NumberOfComponents="3" format="ascii">'
    do e = 1, surf % n_elems
      write(9,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
            IN_5, surf % elem(e) % nx, surf % elem(e) % ny, surf % elem(e) % nz
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'

    ! End of cell data
    write(9,'(a,a)') IN_3, '</CellData>'

    !------------!
    !            !
    !   Footer   !
    !            !
    !------------!
    write(9,'(a,a)') IN_2, '</Piece>'
    write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
    write(9,'(a,a)') IN_0, '</VTKFile>'
    close(9)
  end if

  end subroutine
