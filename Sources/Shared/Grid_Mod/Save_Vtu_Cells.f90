!==============================================================================!
  subroutine Save_Vtu_Cells(Grid, sub, n_nodes_sub, n_cells_sub)
!------------------------------------------------------------------------------!
!   Writes cells in vtu file format                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: sub, n_nodes_sub, n_cells_sub
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)     :: data_size
  integer         :: c, n, s, i_fac, data_offset, cell_offset, fu, s1, s2
  real            :: dist1, dist2
  integer         :: n_conns, n_polyg
  character(SL)   :: name_out
  character(DL*2) :: str1, str2
!==============================================================================!

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) n_conns = n_conns + abs(Grid % cells_n_nodes(c))
  end do

  ! Count face data for polyhedral cells, you will need it later
  n_polyg = 0
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
      if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                  ! add one for number of polyfaces
        do i_fac = 1, Grid % cells_n_faces(c)  ! add all faces and their nodes
          s = Grid % cells_f(i_fac, c)
          n = Grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end if
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File % Set_Name(name_out, processor=sub, extension='.vtu')
  call File % Open_For_Writing_Binary(name_out, fu)

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
  if(Grid % polyhedral) then

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

  end if  ! is Grid polyhedral

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

  ! Number of nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="GridNumberOfNodes"'      //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Wall distance
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridWallDistance"'       //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP  ! prepare for next

  ! Cell volume
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridCellVolume"'         //  &
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
  data_size = int(n_nodes_sub * RP * 3, SP)
  write(fu) data_size
  do n = 1, Grid % n_nodes
    if(Grid % new_n(n) .ne. 0) then
      write(fu) Grid % xn(n), Grid % yn(n), Grid % zn(n)
    end if
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Cells' nodes
  data_size = int(n_conns * IP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then

      ! Tetrahedral, pyramid, wedge and hexahedral cells
      if( any( Grid % cells_n_nodes(c) .eq. (/4,5,6,8/)  ) ) then
        write(fu) Grid % new_n(Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

      ! Polyhedral cells
      else if(Grid % cells_n_nodes(c) < 0) then
        write(fu) Grid % new_n(Grid % cells_n(1:-Grid % cells_n_nodes(c), c))-1

      else
        print *, '# Unsupported cell type with ',  &
                    Grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end if
  end do

  ! Cells' offsets
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  cell_offset = 0
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      cell_offset = cell_offset + abs(Grid % cells_n_nodes(c))
      write(fu) cell_offset
    end if
  end do

  ! Cells' types
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      if(Grid % cells_n_nodes(c) .eq. 4) write(fu) VTK_TETRA
      if(Grid % cells_n_nodes(c) .eq. 8) write(fu) VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c) .eq. 6) write(fu) VTK_WEDGE
      if(Grid % cells_n_nodes(c) .eq. 5) write(fu) VTK_PYRAMID
      if(Grid % cells_n_nodes(c) .lt. 0) write(fu) VTK_POLYHEDRON
    end if
  end do

  ! For polyhedral grids, save faces and face offsets
  if(Grid % polyhedral) then

    ! Write polyhedral cells' faces
    data_size = int(n_polyg * IP, SP)
    write(fu) data_size
    do c = 1, Grid % n_cells
      if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          write(fu) Grid % cells_n_faces(c)      ! write number of its polyfaces
          do i_fac = 1, Grid % cells_n_faces(c)  ! and all polyfaces
            s = Grid % cells_f(i_fac, c)
            if(Grid % faces_s(s) .ne. 0) then    ! face has a shadow, if it ...
              s1 = s                             ! ... is closer, plot that!
              s2 = Grid % faces_s(s)
              dist1 = Math % Distance(                              &
                      Grid % xc(c),  Grid % yc(c),  Grid % zc(c),   &
                      Grid % xf(s1), Grid % yf(s1), Grid % zf(s1))
              dist2 = Math % Distance(                              &
                      Grid % xc(c),  Grid % yc(c),  Grid % zc(c),   &
                      Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
              if(dist1 < dist2) s = s1
              if(dist2 < dist1) s = s2
            end if
            n = Grid % faces_n_nodes(s)
            write(fu) n, Grid % new_n(Grid % faces_n(1:n, s))-1
          end do
        end if
      end if
    end do

    ! Write polyhedral cells' faces offsets
    data_size = int(Grid % n_cells * IP, SP)
    write(fu) data_size
    cell_offset = 0
    do c = 1, Grid % n_cells
      if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          cell_offset = cell_offset + 1          ! to store number of polyfaces
          do i_fac = 1, Grid % cells_n_faces(c)  ! to store polyfaces
            s = Grid % cells_f(i_fac, c)
            n = Grid % faces_n_nodes(s)
            cell_offset = cell_offset + 1 + n    ! number of nodes and nodes
          end do
          write(fu) cell_offset                  ! write the current offset
        else
          write(fu) -1             ! not a polyhedron, offsets are not needed
        end if
      end if
    end do

  end if  ! if Grid % polyhedral

  !---------------!
  !   Cell data   !
  !---------------!

  ! Processor i.d.
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      write(fu) Grid % Comm % cell_proc(c)
    end if
  end do

  ! Number of nodes
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      write(fu) abs(Grid % cells_n_nodes(c))
    end if
  end do

  ! Wall distance
  data_size = int(n_cells_sub * RP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      write(fu) Grid % wall_dist(c)
    end if
  end do

  ! Cell volume
  data_size = int(n_cells_sub * RP, SP)
  write(fu) data_size
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      write(fu) Grid % vol(c)
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
  if(maxval(Grid % Comm % cell_proc(:)) > 1 .and. sub .eq. 1) then

    call File % Set_Name(name_out, extension='.pvtu')
    call File % Open_For_Writing_Ascii(name_out, fu)

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
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" ' //   &
                            ' Name="GridWallDistance"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type="Float64" ' //  &
                            ' Name="GridCellVolume"/>'
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, maxval(Grid % Comm % cell_proc(:))
      call File % Set_Name(name_out, processor=n, extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  end subroutine
