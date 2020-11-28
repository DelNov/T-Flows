!==============================================================================!
  subroutine Save_Vtu_Cells(grid, sub, n_nodes_sub, n_cells_sub)
!------------------------------------------------------------------------------!
!   Writes: name.vtu, name.faces.vtu, name.shadow.vtu                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, n_nodes_sub, n_cells_sub
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)     :: data_size
  integer         :: c, n, s, i_pol, data_offset, cell_offset, fu
  integer         :: n_conns, n_polyg
  character(SL)   :: name_out
  character(DL*2) :: str1, str2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP = DP  ! int. precision is double precision
  integer, parameter :: RP = DP  ! real precision is double precision
!==============================================================================!

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) n_conns = n_conns + abs(grid % cells_n_nodes(c))
  end do

  ! Count face data for polyhedral cells, you will need it later
  n_polyg = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
      if(grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                  ! add one for number of polyfaces
        do i_pol = 1, grid % cells_n_polyg(c)  ! add all faces and their nodes
          s = grid % cells_p(i_pol, c)
          n = grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end if
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File_Mod_Set_Name(name_out, processor=sub, extension='.vtu')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>'             // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid"'  //  &
                    ' version="0.1"'                    //  &
                    ' byte_order="LittleEndian">'       // LF
  write(fu) IN_1 // '<UnstructuredGrid>' // LF
  write(str1, '(i0.0)') n_nodes_sub
  write(str2, '(i0.0)') n_cells_sub
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"' //  &
                    ' NumberOfCells="' // trim(str2) // '">'       // LF
  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + n_nodes_sub * RP * 3  ! prepare for next

  !-----------!
  !   Cells   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Cells' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Cells' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! For polyhedral grids, save faces and face offsets
  if(grid % polyhedral) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next

    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                      ' Name="faceoffsets"'            //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  end if

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !---------------!
  !   Cell data   !
  !---------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  ! Processor i.d.
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="Processor"'              //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Wall distance
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="WallDistance"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP  ! prepare for next

  ! Cell volume
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="CellVolume"'             //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP  ! prepare for next

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu) IN_3 // '</CellData>'         // LF
  write(fu) IN_2 // '</Piece>'            // LF
  write(fu) IN_1 // '</UnstructuredGrid>' // LF

  !-------------------!
  !                   !
  !   Appended data   !
  !                   !
  !-------------------!
  write(fu) IN_0 // '<AppendedData encoding="raw">' // LF
  write(fu) '_'

  !-----------!
  !   Nodes   !
  !-----------!
  data_size = n_nodes_sub * RP * 3
  write(fu) data_size
  do n = 1, grid % n_nodes
    if(grid % new_n(n) .ne. 0)  &
      write(fu) grid % xn(n), grid % yn(n), grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Cells' nodes
  data_size = n_conns * IP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then

      ! Tetrahedral, pyramid, wedge and hexahedral cells
      if( any( grid % cells_n_nodes(c) .eq. (/4,5,6,8/)  ) ) then
        write(fu) grid % new_n(grid % cells_n(1:grid % cells_n_nodes(c), c))-1

      ! Polyhedral cells
      else if(grid % cells_n_nodes(c) < 0) then
        write(fu) grid % new_n(grid % cells_n(1:-grid % cells_n_nodes(c), c))-1

      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end if
  end do

  ! Cells' offsets
  data_size = n_cells_sub * IP
  write(fu) data_size
  cell_offset = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      cell_offset = cell_offset + abs(grid % cells_n_nodes(c))
      write(fu) cell_offset
    end if
  end do

  ! Cells' types
  data_size = n_cells_sub * IP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      if(grid % cells_n_nodes(c) .eq. 4) write(fu) VTK_TETRA
      if(grid % cells_n_nodes(c) .eq. 8) write(fu) VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c) .eq. 6) write(fu) VTK_WEDGE
      if(grid % cells_n_nodes(c) .eq. 5) write(fu) VTK_PYRAMID
      if(grid % cells_n_nodes(c) .lt. 0) write(fu) VTK_POLYHEDRON
    end if
  end do

  ! For polyhedral grids, save faces and face offsets
  if(grid % polyhedral) then

    ! Write polyhedral cells' faces
    data_size = n_polyg * IP
    write(fu) data_size
    do c = 1, grid % n_cells
      if(grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          write(fu) grid % cells_n_polyg(c)      ! write number of its polyfaces
          do i_pol = 1, grid % cells_n_polyg(c)  ! and all polyfaces
            s = grid % cells_p(i_pol, c)
            n = grid % faces_n_nodes(s)
            write(fu) n, grid % new_f(grid % faces_n(1:n, s))-1
          end do
        end if
      end if
    end do

    ! Write polyhedral cells' faces offsets
    data_size = grid % n_cells * IP
    write(fu) data_size
    cell_offset = 0
    do c = 1, grid % n_cells
      if(grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          cell_offset = cell_offset + 1          ! to store number of polyfaces
          do i_pol = 1, grid % cells_n_polyg(c)  ! to store polyfaces
            s = grid % cells_p(i_pol, c)
            n = grid % faces_n_nodes(s)
            cell_offset = cell_offset + 1 + n    ! number of nodes and nodes
          end do
          write(fu) cell_offset                  ! write the current offset
        else
          write(fu) -1             ! not a polyhedron, offsets are not needed
        end if
      end if
    end do

  end if

  !---------------!
  !   Cell data   !
  !---------------!

  ! Processor i.d.
  data_size = n_cells_sub * IP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % comm % cell_proc(c)
    end if
  end do

  ! Wall distance
  data_size = n_cells_sub * RP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % wall_dist(c)
    end if
  end do

  ! Cell volume
  data_size = n_cells_sub * RP
  write(fu) data_size
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % vol(c)
    end if
  end do

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!

  ! Create it only from subdomain 1, when decomposed
  if(maxval(grid % comm % cell_proc(:)) > 1 .and. sub .eq. 1) then

    call File_Mod_Set_Name(name_out, extension='.pvtu')
    call File_Mod_Open_File_For_Writing(name_out, fu)

    ! Header
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(fu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(fu,'(a,a)') IN_2, '<PPoints>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64"'// &
                           ' NumberOfComponents="3"/>'
    write(fu,'(a,a)') IN_2, '</PPoints>'

    ! Data section is not mandatory, but very useful
    write(fu,'(a,a)') IN_2, '<PCellData Scalars="scalars" vectors="velocity">'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Int64" Name="Processor"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" Name="WallDistance"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" Name="CellVolume"/>'
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, maxval(grid % comm % cell_proc(:))
      call File_Mod_Set_Name(name_out, processor=n, extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  end subroutine
