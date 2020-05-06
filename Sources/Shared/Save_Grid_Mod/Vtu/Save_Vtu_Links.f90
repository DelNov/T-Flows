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
  integer           :: n, c, c1, c2, s, offset, fu
  integer           :: nf_sub_non_per, nf_sub_per
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, processor=sub, extension='.links.vtu')
  call File_Mod_Open_File_For_Writing(name_out, fu)

  !-----------!
  !   Start   !
  !-----------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//  &
                         'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                    IN_2, '<Piece NumberOfPoints="',               &
                           n_nodes_sub + n_cells_sub +             &
                           n_bnd_cells_sub + n_buf_cells_sub,      &
                          '" NumberOfCells ="',                    &
                           n_cells_sub + n_faces_sub + n_buf_cells_sub, '">' 
  !-----------!
  !   Nodes   !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                          '="3" format="ascii">'
  do n = 1, grid % n_nodes
    if(grid % new_n(n) .ne. 0) write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')   &
                   IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub)      &
      write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')   &
                 IN_5, grid % xc(c), grid % yc(c), grid % zc(c)
  end do
  do c = -1,-grid % n_bnd_cells,-1
    if(grid % comm % cell_proc(c) .eq. sub)      &
      write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)')   &
                 IN_5, grid % xc(c), grid % yc(c), grid % zc(c)
  end do
  do c = 1,n_buf_cells_sub
    write(fu, '(a,1pe15.7,1pe15.7,1pe15.7)') IN_5,  &
                grid % xc(buf_recv_ind(c)),         &
                grid % yc(buf_recv_ind(c)),         &
                grid % zc(buf_recv_ind(c))
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
    if(grid % comm % cell_proc(c) .eq. sub) then

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
          write(fu,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  & 
                                    n_nodes_sub + c2 - 1
        else
          write(fu,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  &
                                    n_nodes_sub + n_cells_sub - c2 - 1
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
          write(fu,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  &
                                    n_nodes_sub + c2 - 1
        else
          write(fu,'(a,2i9)') IN_5, n_nodes_sub + c1 - 1,  &
                                    n_nodes_sub + n_cells_sub - c2 - 1
        end if
      end if

    end if
  end do  

  ! Interprocessor links
  do c = 1, n_buf_cells_sub
    c1 = buf_send_ind(c) 
    write(fu,'(a,2i9)') IN_5,  &
      n_nodes_sub + c1 - 1, n_nodes_sub + n_cells_sub + n_bnd_cells_sub + c - 1
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'

  print '(a38,i9)', '# Non-periodic links    :            ', nf_sub_non_per
  print '(a38,i9)', '# Periodic links        :            ', nf_sub_per
  print '(a38,i9)', '# Inter-processor links :            ', n_buf_cells_sub

  ! Now write all cells' offsets
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets" ' // & 
                          'format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      offset = offset + grid % cells_n_nodes(c)
      write(fu,'(a,i9)') IN_5, offset
    end if
  end do
  do c = 1, nf_sub_non_per
    offset = offset + 2
    write(fu,'(a,i9)') IN_5, offset
  end do
  do c = 1, nf_sub_per
    offset = offset + 2
    write(fu,'(a,i9)') IN_5, offset
  end do
  do c = 1, n_buf_cells_sub
    offset = offset + 2
    write(fu,'(a,i9)') IN_5, offset
  end do

  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      if(grid % cells_n_nodes(c) .eq. 4) write(fu,'(a,i9)') IN_5, VTK_TETRA
      if(grid % cells_n_nodes(c) .eq. 8) write(fu,'(a,i9)') IN_5, VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c) .eq. 6) write(fu,'(a,i9)') IN_5, VTK_WEDGE
      if(grid % cells_n_nodes(c) .eq. 5) write(fu,'(a,i9)') IN_5, VTK_PYRAMID
    end if
  end do
  do c = 1, nf_sub_non_per
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do
  do c = 1, nf_sub_per
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do
  do c = 1, n_buf_cells_sub
    write(fu,'(a,i9)') IN_5, VTK_LINE
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !----------------!
  !   Link types   !
  !----------------!
  write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'
  write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' // & 
                          'Name="link type" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu,'(a,i9)') IN_5, 0
    end if
  end do
  do c = 1, nf_sub_non_per
    write(fu,'(a,i9)') IN_5, 1
  end do
  do c = 1, nf_sub_per
    write(fu,'(a,i9)') IN_5, 2
  end do
  do c = 1, n_buf_cells_sub
    write(fu,'(a,i9)') IN_5, 3
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</CellData>'

  !------------!
  !   Footer   !
  !------------!
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'

  close(fu)

  end subroutine
