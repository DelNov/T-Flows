!==============================================================================!
  subroutine Save_Vtu_Links(grid,             &
                            sub,              &  ! subdomain
                            n_nodes_sub,      &  ! number of nodes in the sub. 
                            n_cells_sub,      &  ! number of cells in the sub. 
                            n_faces_sub,      &  ! number of faces in the sub.
                            n_bnd_cells_sub,  &  ! number of bnd. cells in sub
                            n_buf_cells_sub)     ! number of buf. cells in sub
!------------------------------------------------------------------------------!
!   Creates the file "name.ln.vtu" to check the cell connections.              !
!                                                                              !
!   Links between the computational cells have been introduced as aditional    !
!   cells of general type. Cell centers are introduced as aditional nodes.     !
!   Material of these links is different than from the cells, so that they     !
!   can be visualised  more easily in GMV.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, n_nodes_sub, n_cells_sub, n_faces_sub,  &
                          n_bnd_cells_sub, n_buf_cells_sub
!-----------------------------------[Locals]-----------------------------------!
  integer(4)         :: data_size
  integer            :: n, c, c1, c2, s, data_offset, offset, fu
  integer            :: nf_sub_non_per, nf_sub_per
  integer            :: n_nodes_here, n_cells_here
  character(len=160) :: name_out, str1, str2
  integer, parameter :: IP=8, RP=8, SP=4
!==============================================================================!

  n_nodes_here = n_nodes_sub + n_cells_sub +             &
                 n_bnd_cells_sub + n_buf_cells_sub
  n_cells_here = n_cells_sub + n_faces_sub + n_buf_cells_sub

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, processor=sub, extension='.links.vtu')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>' // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid" version="0.1" '//  &
                    ' byte_order="LittleEndian">' // LF
  write(fu) IN_1 // '<UnstructuredGrid>' // LF
  write(str1, '(i0.0)') n_nodes_here
  write(str2, '(i0.0)') n_cells_here
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
  data_offset = data_offset + SP + n_nodes_here * RP * 3  ! prepare for next

  !-----------!
  !   Cells   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! First write all cells' nodes
  write(fu) IN_4 // '<DataArray type="Int64" Name="connectivity"' //  &
                    ' format="ascii">' // LF

  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then

      ! Hexahedral
      if(grid % cells_n_nodes(c) .eq. 8) then
        write(str1,'(a,8i9)')                      &
           IN_5,                                 &
           grid % new_n(grid % cells_n(1,c))-1,  &
           grid % new_n(grid % cells_n(2,c))-1,  &
           grid % new_n(grid % cells_n(4,c))-1,  &
           grid % new_n(grid % cells_n(3,c))-1,  &
           grid % new_n(grid % cells_n(5,c))-1,  &
           grid % new_n(grid % cells_n(6,c))-1,  &
           grid % new_n(grid % cells_n(8,c))-1,  &
           grid % new_n(grid % cells_n(7,c))-1
        write(fu) trim(str1) // LF

      ! Wedge
      else if(grid % cells_n_nodes(c) .eq. 6) then
        write(str1,'(a,6i9)')                      &
           IN_5,                                 &
           grid % new_n(grid % cells_n(1,c))-1,  &
           grid % new_n(grid % cells_n(2,c))-1,  &
           grid % new_n(grid % cells_n(3,c))-1,  &
           grid % new_n(grid % cells_n(4,c))-1,  &
           grid % new_n(grid % cells_n(5,c))-1,  &
           grid % new_n(grid % cells_n(6,c))-1
        write(fu) trim(str1) // LF

      ! Tetrahedra
      else if(grid % cells_n_nodes(c) .eq. 4) then
        write(str1,'(a,4i9)')                      &
           IN_5,                                 &
           grid % new_n(grid % cells_n(1,c))-1,  &
           grid % new_n(grid % cells_n(2,c))-1,  &
           grid % new_n(grid % cells_n(3,c))-1,  &
           grid % new_n(grid % cells_n(4,c))-1
        write(fu) trim(str1) // LF

      ! Pyramid
      else if(grid % cells_n_nodes(c) .eq. 5) then
        write(str1,'(a,5i9)')                      &
           IN_5,                                 &
           grid % new_n(grid % cells_n(1,c))-1,  &
           grid % new_n(grid % cells_n(2,c))-1,  &
           grid % new_n(grid % cells_n(4,c))-1,  &
           grid % new_n(grid % cells_n(3,c))-1,  &
           grid % new_n(grid % cells_n(5,c))-1
        write(fu) trim(str1) // LF

      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end if
  end do

  ! Physical links; non-periodic
  nf_sub_non_per = 0 
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if( (grid % new_f(s) > 0) .and. (grid % new_f(s) <= n_faces_sub) ) then

      if( (grid % sx(s) * (grid % xc(c2)-grid % xc(c1) ) +  &
           grid % sy(s) * (grid % yc(c2)-grid % yc(c1) ) +  &
           grid % sz(s) * (grid % zc(c2)-grid % zc(c1) ))  > 0.0 ) then

        nf_sub_non_per = nf_sub_non_per + 1

        c1 = grid % new_c(grid % faces_c(1,s))
        c2 = grid % new_c(grid % faces_c(2,s))
        if( c2  > 0 ) then
          write(str1,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  & 
                                      n_nodes_sub + c2 - 1
          write(fu) trim(str1) // LF
        else
          write(str1,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  &
                                      n_nodes_sub + n_cells_sub - c2 - 1
          write(fu) trim(str1) // LF
        end if
      end if

    end if
  end do

  ! Physical links; periodic
  nf_sub_per    = 0 
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if( (grid % new_f(s) > 0) .and. (grid % new_f(s) <= n_faces_sub) ) then

      if( (grid % sx(s) * (grid % xc(c2)-grid % xc(c1) ) +  &
           grid % sy(s) * (grid % yc(c2)-grid % yc(c1) ) +  &
           grid % sz(s) * (grid % zc(c2)-grid % zc(c1) ))  < 0.0 ) then

        nf_sub_per = nf_sub_per + 1

        c1 = grid % new_c(grid % faces_c(1,s))
        c2 = grid % new_c(grid % faces_c(2,s))
        if( c2  > 0 ) then
          write(str1,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  & 
                                      n_nodes_sub + c2 - 1
          write(fu) trim(str1) // LF
        else
          write(str1,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  &
                                      n_nodes_sub + n_cells_sub - c2 - 1
          write(fu) trim(str1) // LF
        end if
      end if

    end if
  end do  

  ! Interprocessor links
  do c = 1, n_buf_cells_sub
    c1 = buf_send_ind(c) 
    write(str1,'(a,2i9)') IN_5,  &
      n_nodes_sub + c1 - 1, n_nodes_sub + n_cells_sub + n_bnd_cells_sub + c - 1
    write(fu) trim(str1) // LF
  end do

  write(fu) IN_4 // '</DataArray>' // LF

  print '(a38,i9)', '# Non-periodic links    :            ', nf_sub_non_per
  print '(a38,i9)', '# Periodic links        :            ', nf_sub_per
  print '(a38,i9)', '# Inter-processor links :            ', n_buf_cells_sub

  ! Now write all cells' offsets
  write(fu) IN_4 // '<DataArray type="Int64" Name="offsets" ' // & 
                    'format="ascii">' // LF
  offset = 0
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      offset = offset + grid % cells_n_nodes(c)
      write(str1,'(a,i9)') IN_5, offset
      write(fu) trim(str1) // LF
    end if
  end do
  do c = 1, nf_sub_non_per
    offset = offset + 2
    write(str1,'(a,i9)') IN_5, offset
    write(fu) trim(str1) // LF
  end do
  do c = 1, nf_sub_per
    offset = offset + 2
    write(str1,'(a,i9)') IN_5, offset
    write(fu) trim(str1) // LF
  end do
  do c = 1, n_buf_cells_sub
    offset = offset + 2
    write(str1,'(a,i9)') IN_5, offset
    write(fu) trim(str1) // LF
  end do

  write(fu) IN_4 // '</DataArray>' // LF

  ! Now write all cells' types
  write(fu) IN_4 // '<DataArray type="Int64" Name="types" format="ascii">' // LF
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      if(grid % cells_n_nodes(c) .eq. 4) then
        write(str1,'(a,i9)') IN_5, VTK_TETRA
        write(fu) trim(str1) // LF
      end if
      if(grid % cells_n_nodes(c) .eq. 8) then
        write(str1,'(a,i9)') IN_5, VTK_HEXAHEDRON
        write(fu) trim(str1) // LF
      end if
      if(grid % cells_n_nodes(c) .eq. 6) then
        write(str1,'(a,i9)') IN_5, VTK_WEDGE
        write(fu) trim(str1) // LF
      end if
      if(grid % cells_n_nodes(c) .eq. 5) then
        write(str1,'(a,i9)') IN_5, VTK_PYRAMID
        write(fu) trim(str1) // LF
      end if
    end if
  end do
  do c = 1, nf_sub_non_per
    write(str1,'(a,i9)') IN_5, VTK_LINE
    write(fu) trim(str1) // LF
  end do
  do c = 1, nf_sub_per
    write(str1,'(a,i9)') IN_5, VTK_LINE
    write(fu) trim(str1) // LF
  end do
  do c = 1, n_buf_cells_sub
    write(str1,'(a,i9)') IN_5, VTK_LINE
    write(fu) trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Cells>'     // LF

  !----------------!
  !   Link types   !
  !----------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF
  write(fu) IN_4 // '<DataArray type="Int64" ' // & 
                          'Name="link type" format="ascii">' // LF
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(str1,'(a,i9)') IN_5, 0
      write(fu) trim(str1) // LF
    end if
  end do
  do c = 1, nf_sub_non_per
    write(str1,'(a,i9)') IN_5, 1
    write(fu) trim(str1) // LF
  end do
  do c = 1, nf_sub_per
    write(str1,'(a,i9)') IN_5, 2
    write(fu) trim(str1) // LF
  end do
  do c = 1, n_buf_cells_sub
    write(str1,'(a,i9)') IN_5, 3
    write(fu) trim(str1) // LF
  end do
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</CellData>'  // LF

  !------------!
  !   Footer   !
  !------------!
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
  data_size = n_nodes_here * RP * 3
  write(fu) data_size
  do n = 1, grid % n_nodes
    if(grid % new_n(n) .ne. 0) then
      write(fu) grid % xn(n), grid % yn(n), grid % zn(n)
    end if
  end do
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % xc(c), grid % yc(c), grid % zc(c)
    end if
  end do
  do c = -1,-grid % n_bnd_cells,-1
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % xc(c), grid % yc(c), grid % zc(c)
    end if
  end do
  do c = 1,n_buf_cells_sub
    write(fu) grid % xc(buf_recv_ind(c)),  &
              grid % yc(buf_recv_ind(c)),  &
              grid % zc(buf_recv_ind(c))
  end do

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>'            // LF

  close(fu)

  end subroutine
