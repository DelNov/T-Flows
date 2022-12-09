!==============================================================================!
  subroutine Save_Cfn(Grid,        &
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
  class(Grid_Type)    :: Grid
  integer, intent(in) :: sub, nn_sub, nc_sub, nf_sub, ns_sub, nbc_sub
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, fu, c1, c2, ss, sr
  character(SL) :: name_out
!==============================================================================!

  call Profiler % Start('Save_Cfn')

  !----------------------!
  !   Create .cfn file   !
  !----------------------!
  call File % Set_Name(name_out,         &
                       processor=sub,    &
                       extension='.cfn')
  call File % Open_For_Writing_Binary(name_out, fu)

  !-------------------------!
  !   Save real precision   !
  !-------------------------!
  write(fu) RP

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  write(fu) nn_sub
  write(fu) nc_sub              ! new way: add buffer cells to cells
  write(fu) nbc_sub             ! number of boundary cells
  write(fu) nf_sub
  write(fu) ns_sub
  write(fu) Grid % n_bnd_cond  ! number of bounary conditions

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  write(fu) Grid % polyhedral

  !------------------------!
  !   Domain (grid) name   !
  !------------------------!
  write(fu) Grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, Grid % n_bnd_cond
    write(fu) Grid % bnd_cond % name(n)
  end do

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  do n = 1, Grid % n_nodes
    if(Grid % new_n(n) > 0) then
      write(fu) Grid % comm % node_glo(n)
    end if
  end do

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) Grid % cells_n_nodes(Grid % old_c(c))
    end if
  end do

  ! Error trap for number of nodes for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(Grid % cells_n_nodes(Grid % old_c(c)) .eq. 0) then
          print *, '# ERROR: Number of nodes is zero at cell:', Grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Cells' nodes
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do n = 1, abs(Grid % cells_n_nodes(Grid % old_c(c)))
        write(fu) Grid % new_n(Grid % cells_n(n, Grid % old_c(c)))
      end do
    end if
  end do

  ! Error trap for cells' nodes
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do n = 1, abs(Grid % cells_n_nodes(Grid % old_c(c)))
        if(Grid % new_n(Grid % cells_n(n, Grid % old_c(c))) .eq. 0) then
          print *, '# ERROR: Node index is zero at cell:', Grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Number of faces for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) Grid % cells_n_faces(Grid % old_c(c))
    end if
  end do

  ! Error trap for number of faces for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(Grid % cells_n_faces(Grid % old_c(c)) .eq. 0) then
          print *, '# ERROR: Number of faces is zero at cell:', Grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Cells' faces.  They are still faces, kept in the same array.
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do s = 1, Grid % cells_n_faces(Grid % old_c(c))
        write(fu) Grid % new_f(Grid % cells_f(s, Grid % old_c(c)))
      end do
    end if
  end do

  ! Error trap for cells' faces
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do s = 1, Grid % cells_n_faces(Grid % old_c(c))
        if(Grid % new_f(Grid % cells_f(s, Grid % old_c(c))) .eq. 0) then
          print *, '# ERROR: Face index is zero at cell:', Grid % old_c(c)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Cells' processor ids
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) Grid % comm % cell_proc(Grid % old_c(c))
    end if
  end do

  ! Cells' global indices
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) Grid % comm % cell_glo(Grid % old_c(c))
    end if
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      write(fu) Grid % faces_n_nodes(Grid % old_f(s))
    end if
  end do

  ! Error trap for number of nodes for each face
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      if(Grid % faces_n_nodes(Grid % old_f(s)) .eq. 0) then
        print *, '# ERROR: Number of nodes is zero at face:', Grid % old_f(s)
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  ! Faces' nodes
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      do n = 1, Grid % faces_n_nodes(Grid % old_f(s))
        write(fu) Grid % new_n(Grid % faces_n(n, Grid % old_f(s)))
      end do
    end if
  end do

  ! Error trap for faces' nodes
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      do n = 1, Grid % faces_n_nodes(Grid % old_f(s))
        if(Grid % new_n(Grid % faces_n(n, Grid % old_f(s))) .eq. 0) then
          print *, '# ERROR: Node index is zero at face:', Grid % old_f(s)
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end do
    end if
  end do

  ! Faces' cells
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      c1 = Grid % faces_c(1, Grid % old_f(s))
      c2 = Grid % faces_c(2, Grid % old_f(s))

      ! At least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(Grid % comm % cell_proc(c1) .eq. sub .or.  &
         Grid % comm % cell_proc(c2) .eq. sub) then
        if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
          write(fu) Grid % new_c(c1), Grid % new_c(c2)
        else
          write(fu) Grid % new_c(c2), Grid % new_c(c1)
        end if

      ! Face is not active
      else
        write(fu) 0, 0
      end if
    end if
  end do

  ! Error trap for faces' cells
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      c1 = Grid % faces_c(1, Grid % old_f(s))
      c2 = Grid % faces_c(2, Grid % old_f(s))

      ! Check only if least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(Grid % comm % cell_proc(c1) .eq. sub .or.  &
         Grid % comm % cell_proc(c2) .eq. sub) then
        if(Grid % new_c(c1) .eq. 0) then
          print *, '# ERROR: Cell one is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
        if(Grid % new_c(c2) .eq. 0) then
          print *, '# ERROR: Cell two is zero at face:', s
          print *, '# This error is critical.  Exiting!'
          call Comm_Mod_End
          stop
        end if
      end if
    end if
  end do

  ! Faces' shadows
  do sr = 1, Grid % n_faces
    if(Grid % old_f(sr) .ne. 0) then         ! this face will be saved
      ss = Grid % faces_s(Grid % old_f(sr))  ! fetch its shadow

      if(ss .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        write(fu) Grid % new_f(ss)
      else
        write(fu) 0
      end if
    end if
  end do

  do ss = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(ss) .ne. 0) then  ! this shadow face will be saved
      sr = Grid % faces_s(ss)         ! fetch its real counterpart

      if(sr .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        write(fu) Grid % new_f(sr)
      else
        print *, '# ERROR: Shadow faces points to zero face'
        print *, '# This error is critical.  Exiting!'
        call Comm_Mod_End
        stop
      end if
    end if
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  do c = -Grid % n_bnd_cells, -1
    if(Grid % old_c(c) .ne. 0) then
      write(fu) Grid % bnd_cond % color(Grid % old_c(c))
    end if
  end do

  close(fu)

  call Profiler % Stop('Save_Cfn')

  end subroutine
