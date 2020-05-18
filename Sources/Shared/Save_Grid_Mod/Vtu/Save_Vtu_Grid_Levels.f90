!==============================================================================!
  subroutine Save_Vtu_Grid_Levels(grid, mg_lev)
!------------------------------------------------------------------------------!
!   Writes .mg00.faces.vtu file.                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: mg_lev
!-----------------------------------[Locals]-----------------------------------!
  integer            :: c, c1, c2, n, s, offset, level_n_faces, fu
  character(len=80)  :: extension, name_out
!==============================================================================!

  !------------------------------------------------------!
  !   Count number of faces for plotting on this level   !
  !------------------------------------------------------!
  level_n_faces = 0
  do s = 1, grid % n_faces
    call Grid_Mod_Get_C1_And_C2_At_Level(grid, mg_lev, s, c1, c2)
    if(c1 .ne. c2 .and. c2 > 0) level_n_faces = level_n_faces + 1
  end do

  !----------------------!
  !   Create .vtu file   !
  !----------------------!
  extension = '.mg00.faces.vtu'
  write(extension(4:5), '(i2.2)') mg_lev
  call File_Mod_Set_Name(name_out, extension=trim(extension))
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
                    IN_2, '<Piece NumberOfPoints="',       &
                          grid % n_nodes                   &
                        + grid % level(mg_lev) % n_cells,  &
                          '" NumberOfCells ="',            &
                          level_n_faces                    &
                        + grid % level(mg_lev) % n_faces,  &
                          '">'
  !------------!
  !            !
  !   Points   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents=' //  &
                          '"3" format="ascii">'

  !----------------------------------------!
  !   Physical nodes on the finest level   !
  !----------------------------------------!
  do n = 1, grid % n_nodes
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')                &
               IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do

  !---------------------------------------!
  !   Cell centers on the current level   !
  !---------------------------------------!
  do c = 1, grid % level(mg_lev) % n_cells
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')         &
               IN_5, grid % level(mg_lev) % xc(c),  &
                     grid % level(mg_lev) % yc(c),  &
                     grid % level(mg_lev) % zc(c)
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !------------------!
  !                  !
  !   Connectivity   !
  !                  !
  !------------------!
  write(fu,'(a,a)') IN_3, '<Cells>'

  !----------------------------------!
  !   First write all faces' nodes   !
  !----------------------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                          ' format="ascii">'
  do s = 1, grid % n_faces

    call Grid_Mod_Get_C1_And_C2_At_Level(grid, mg_lev, s, c1, c2)

    ! Either boundary on the finest mg_lev,
    ! ... or inter-cell face on any mg_lev.
    if(c1 .ne. c2 .and. c2 > 0) then
      if(grid % faces_n_nodes(s) .eq. 4) then
        write(fu,'(a,4i9)')                               &
          IN_5,                                          &
          grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
          grid % faces_n(3,s)-1, grid % faces_n(4,s)-1
      else if(grid % faces_n_nodes(s) .eq. 3) then
        write(fu,'(a,3I9)')                               &
          IN_5,                                          &
          grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
          grid % faces_n(3,s)-1
      else
        print *, '# Unsupported cell type ',       &
                   grid % faces_n_nodes(s), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end if
  end do

  !-------------------------------!
  !   The face-cell connections   !
  !-------------------------------!
  do s = 1, grid % level(mg_lev) % n_faces
    write(fu,'(a,2i9)') IN_5,   grid % level(mg_lev) % faces_c(1,s)  &
                              + grid % n_nodes - 1,                  &
                                grid % level(mg_lev) % faces_c(2,s)  &
                              + grid % n_nodes - 1
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'

  !-------------!
  !             !
  !   Offsets   !
  !             !
  !-------------!
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets" ' // &
                          'format="ascii">'
  offset = 0

  !----------------------------!
  !   Offsets for cell faces   !
  !----------------------------!
  do s = 1, grid % n_faces

    call Grid_Mod_Get_C1_And_C2_At_Level(grid, mg_lev, s, c1, c2)

    ! Either boundary on the finest mg_lev,
    ! ... or inter-cell face on any mg_lev.
    if(c1 .ne. c2 .and. c2 > 0) then
      offset = offset + grid % faces_n_nodes(s)
      write(fu,'(a,i9)') IN_5, offset
    end if

  end do

  !----------------------------------!
  !   Face-cell connection offsets   !
  !----------------------------------!
  do s = 1, grid % level(mg_lev) % n_faces
    offset = offset + 2
    write(fu,'(a,i9)') IN_5, offset
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'

  !------------------!
  !                  !
  !   Cells' types   !
  !                  !
  !------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'

  !-----------!
  !   Faces   !
  !-----------!
  do s = 1, grid % n_faces

    call Grid_Mod_Get_C1_And_C2_At_Level(grid, mg_lev, s, c1, c2)

    ! Either boundary on the finest mg_lev,
    ! ... or inter-cell face on any mg_lev.
    if(c1 .ne. c2 .and. c2 > 0) then
      if(grid % faces_n_nodes(s) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_QUAD
      if(grid % faces_n_nodes(s) .eq. 3) write(fu,'(a,i9)') IN_5, VTK_TRIANGLE
    end if

  end do

  !-----------------!
  !   Connections   !
  !-----------------!
  do s = 1, grid % level(mg_lev) % n_faces
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !               !
  !   Face data   !
  !               !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  ! Boundary conditions
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" '   // &
                   'Name="CellType" format="ascii">'
  do s = 1, grid % n_faces

    call Grid_Mod_Get_C1_And_C2_At_Level(grid, mg_lev, s, c1, c2)

    if(c1 .ne. c2 .and. c2 > 0) then
      write(fu,'(a,i9)') IN_5, 0
    end if
  end do

  do s = 1, grid % level(mg_lev) % n_faces
    write(fu,'(a,i9)') IN_5, -1
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_3, '</CellData>'
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'

  close(fu)

  end subroutine
