!==============================================================================!
  subroutine Grid_Mod_Save_Cfn(grid,        &
                               sub,         &  ! subdomain
                               nn_sub,      &  ! number of nodes in the sub. 
                               nc_sub,      &  ! number of cells in the sub. 
                               nf_sub,      &  ! number of faces in the sub.
                               ns_sub,      &  ! number of shadow faces
                               nbc_sub)
!------------------------------------------------------------------------------!
!   Writes file with connectivity data: name.cfn                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, nn_sub, nc_sub, nf_sub, ns_sub, nbc_sub
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, fu, c1, c2, ss, sr
  character(SL) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .cfn file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out,         &
                         processor=sub,    &
                         extension='.cfn')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  write(fu) nn_sub
  write(fu) nc_sub              ! new way: add buffer cells to cells
  write(fu) nbc_sub             ! number of boundary cells
  write(fu) nf_sub
  write(fu) ns_sub
  write(fu) grid % n_bnd_cond  ! number of bounary conditions

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  write(fu) grid % polyhedral

  !------------------------!
  !   Domain (grid) name   !
  !------------------------!
  write(fu) grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, grid % n_bnd_cond
    write(fu) grid % bnd_cond % name(n)
  end do

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  do n = 1, grid % n_nodes
    if(grid % new_n(n) > 0) then
      write(fu) grid % comm % node_glo(n)
    end if
  end do

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) grid % cells_n_nodes(grid % old_c(c))
    end if
  end do

  ! Error trap for number of nodes for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(grid % cells_n_nodes(grid % old_c(c)) .eq. 0) then
          print *, '# ERROR: Number of nodes is zero at cell:', grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Cells' nodes
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do n = 1, abs(grid % cells_n_nodes(grid % old_c(c)))
        write(fu) grid % new_n(grid % cells_n(n, grid % old_c(c)))
      end do
    end if
  end do

  ! Error trap for cells' nodes
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do n = 1, abs(grid % cells_n_nodes(grid % old_c(c)))
        if(grid % new_n(grid % cells_n(n, grid % old_c(c))) .eq. 0) then
          print *, '# ERROR: Node index is zero at cell:', grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Number of faces for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) grid % cells_n_faces(grid % old_c(c))
    end if
  end do

  ! Error trap for number of faces for each cell
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(grid % cells_n_faces(grid % old_c(c)) .eq. 0) then
          print *, '# ERROR: Number of faces is zero at cell:', grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Cells' faces.  They are still faces, kept in the same array.
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do s = 1, grid % cells_n_faces(grid % old_c(c))
        write(fu) grid % new_f(grid % cells_f(s, grid % old_c(c)))
      end do
    end if
  end do

  ! Error trap for cells' faces
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do s = 1, grid % cells_n_faces(grid % old_c(c))
        if(grid % new_f(grid % cells_f(s, grid % old_c(c))) .eq. 0) then
          print *, '# ERROR: Face index is zero at cell:', grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Cells' processor ids
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) grid % comm % cell_proc(grid % old_c(c))
    end if
  end do

  ! Cells' global indices
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) grid % comm % cell_glo(grid % old_c(c))
    end if
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      write(fu) grid % faces_n_nodes(grid % old_f(s))
    end if
  end do

  ! Error trap for number of nodes for each face
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      if(grid % faces_n_nodes(grid % old_f(s)) .eq. 0) then
        print *, '# ERROR: Number of nodes is zero at face:', grid % old_f(s)
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Faces' nodes
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      do n = 1, grid % faces_n_nodes(grid % old_f(s))
        write(fu) grid % new_n(grid % faces_n(n, grid % old_f(s)))
      end do
    end if
  end do

  ! Error trap for faces' nodes
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      do n = 1, grid % faces_n_nodes(grid % old_f(s))
        if(grid % new_n(grid % faces_n(n, grid % old_f(s))) .eq. 0) then
          print *, '# ERROR: Node index is zero at face:', grid % old_f(s)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Faces' cells
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      c1 = grid % faces_c(1, grid % old_f(s))
      c2 = grid % faces_c(2, grid % old_f(s))

      ! At least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(grid % comm % cell_proc(c1) .eq. sub .or.  &
         grid % comm % cell_proc(c2) .eq. sub) then
        if(grid % new_c(c2) < 0 .or. grid % new_c(c1) < grid % new_c(c2)) then
          write(fu) grid % new_c(c1), grid % new_c(c2)
        else
          write(fu) grid % new_c(c2), grid % new_c(c1)
        end if

      ! Face is not active
      else
        write(fu) 0, 0
      end if
    end if
  end do

  ! Error trap for faces' cells
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      c1 = grid % faces_c(1, grid % old_f(s))
      c2 = grid % faces_c(2, grid % old_f(s))

      ! Check only if least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(grid % comm % cell_proc(c1) .eq. sub .or.  &
         grid % comm % cell_proc(c2) .eq. sub) then
        if(grid % new_c(c1) .eq. 0) then
          print *, '# ERROR: Cell one is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
        if(grid % new_c(c2) .eq. 0) then
          print *, '# ERROR: Cell two is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Faces' shadows
  do sr = 1, grid % n_faces
    if(grid % old_f(sr) .ne. 0) then         ! this face will be saved
      ss = grid % faces_s(grid % old_f(sr))  ! fetch its shadow

      if(ss .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        write(fu) grid % new_f(ss)
      else
        write(fu) 0
      end if
    end if
  end do

  do ss = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(ss) .ne. 0) then  ! this shadow face will be saved
      sr = grid % faces_s(ss)         ! fetch its real counterpart

      if(sr .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        write(fu) grid % new_f(sr)
      else
        print *, '# ERROR: Shadow faces points to zero face'
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Faces' global indices
  do s = 1, grid % n_faces + grid % n_shadows
    if(grid % old_f(s) .ne. 0) then
      write(fu) grid % comm % face_glo(grid % old_f(s))
    end if
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  do c = -grid % n_bnd_cells, -1
    if(grid % old_c(c) .ne. 0) then
      write(fu) grid % bnd_cond % color(grid % old_c(c))
    end if
  end do

  close(fu)

  end subroutine
