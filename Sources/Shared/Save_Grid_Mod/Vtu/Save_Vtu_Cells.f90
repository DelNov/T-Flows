!==============================================================================!
  subroutine Save_Vtu_Cells(grid, sub, n_nodes_sub, n_cells_sub)
!------------------------------------------------------------------------------!
!   Writes: name.vtu, name.faces.vtu, name.shadow.vtu                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, n_nodes_sub, n_cells_sub, lev
!-----------------------------------[Locals]-----------------------------------!
  integer            :: c, n, offset, fu
  character(len=80)  :: name_out
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: VTK_TETRA      = 10  ! cells in VTK format
  integer,           parameter :: VTK_HEXAHEDRON = 12  
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
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
  call File_Mod_Set_Name(name_out, processor=sub, extension='.vtu')
  call File_Mod_Open_File_For_Writing(name_out, fu)

  !-----------!
  !   Start   !
  !-----------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//  &
                          'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                    IN_2, '<Piece NumberOfPoints="', n_nodes_sub,      &
                               '" NumberOfCells ="', n_cells_sub, '">'

  !-----------!
  !   Nodes   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                          '="3" format="ascii">'
  do n = 1, grid % n_nodes
    if(grid % new_n(n) .ne. 0) write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')  &
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
    if(grid % new_c(c) .ne. 0) then

      ! Hexahedral
      if(grid % cells_n_nodes(c) .eq. 8) then
        write(fu,'(a,8i9)')                      &
          IN_5,                                 &
          grid % new_n(grid % cells_n(1,c))-1,  &
          grid % new_n(grid % cells_n(2,c))-1,  &
          grid % new_n(grid % cells_n(4,c))-1,  &
          grid % new_n(grid % cells_n(3,c))-1,  &
          grid % new_n(grid % cells_n(5,c))-1,  &
          grid % new_n(grid % cells_n(6,c))-1,  &
          grid % new_n(grid % cells_n(8,c))-1,  &
          grid % new_n(grid % cells_n(7,c))-1

      ! Wedge       
      else if(grid % cells_n_nodes(c) .eq. 6) then
        write(fu,'(a,6i9)')                      &
          IN_5,                                 &
          grid % new_n(grid % cells_n(1,c))-1,  &
          grid % new_n(grid % cells_n(2,c))-1,  &
          grid % new_n(grid % cells_n(3,c))-1,  &
          grid % new_n(grid % cells_n(4,c))-1,  &
          grid % new_n(grid % cells_n(5,c))-1,  &
          grid % new_n(grid % cells_n(6,c))-1

      ! Tetrahedra  
      else if(grid % cells_n_nodes(c) .eq. 4) then
        write(fu,'(a,4i9)')                      &
          IN_5,                                 &
          grid % new_n(grid % cells_n(1,c))-1,  &
          grid % new_n(grid % cells_n(2,c))-1,  &
          grid % new_n(grid % cells_n(3,c))-1,  &
          grid % new_n(grid % cells_n(4,c))-1

      ! Pyramid     
      else if(grid % cells_n_nodes(c) .eq. 5) then
        write(fu,'(a,5i9)')                      &
          IN_5,                                 &
          grid % new_n(grid % cells_n(1,c))-1,  &
          grid % new_n(grid % cells_n(2,c))-1,  &
          grid % new_n(grid % cells_n(4,c))-1,  &
          grid % new_n(grid % cells_n(3,c))-1,  &
          grid % new_n(grid % cells_n(5,c))-1
      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop 
      end if
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                          'Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      offset = offset + grid % cells_n_nodes(c)
      write(fu,'(a,i9)') IN_5, offset
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      if(grid % cells_n_nodes(c) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_TETRA
      if(grid % cells_n_nodes(c) .eq. 8) write(fu,'(a,i9)') IN_5, VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c) .eq. 6) write(fu,'(a,i9)') IN_5, VTK_WEDGE
      if(grid % cells_n_nodes(c) .eq. 5) write(fu,'(a,i9)') IN_5, VTK_PYRAMID
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  ! Processor i.d.
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' //  &
                          'Name="Processor" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu,'(a,i9)') IN_5, grid % comm % proces(c)
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Coarser grid levels
  do lev = 1, grid % n_levels
    write(fu,'(a,a,i2.2,a)') IN_4, '<DataArray type="Int64" ' // &
                                   'Name="GridLevel', lev, '" format="ascii">'
    do c = 1, grid % n_cells
      if(grid % new_c(c) .ne. 0) then
        write(fu,'(a,i9)') IN_5, grid % level(lev) % cell(c)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
  end do

  ! Wall distance
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
                          'Name="WallDistance" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu,'(a,1pe15.7)') IN_5, grid % wall_dist(c)
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Cell volume
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" ' //  &
                          'Name="CellVolume" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu,'(a,1pe15.7)') IN_5, grid % vol(c)
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

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!

  ! Create it only from subdomain 1, when decomposed
  if(maxval(grid % comm % proces(:)) > 1 .and. sub .eq. 1) then

    call File_Mod_Set_Name(name_out, extension='.pvtu')
    call File_Mod_Open_File_For_Writing(name_out, fu)

    ! Header
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(fu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(fu,'(a,a)') IN_2, '<PPoints>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" NumberOfComponents='// &
                           '"3" format="ascii"/>'
    write(fu,'(a,a)') IN_2, '</PPoints>'

    ! Data section is not mandatory, but very useful
    write(fu,'(a,a)') IN_2, '<PCellData Scalars="scalars" vectors="velocity">'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Int64" Name="Processor"' // &
                           ' format="ascii"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" Name="WallDistance"'// &
                           ' format="ascii"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" Name="CellVolume"' // &
                           ' format="ascii"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" Name="CellDelta"' // &
                           ' format="ascii"/>'
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, maxval(grid % comm % proces(:))
      call File_Mod_Set_Name(name_out, processor=n, extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  end subroutine
