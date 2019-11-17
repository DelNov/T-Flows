!==============================================================================!
  subroutine Grid_Mod_Save_Cns(grid,        &
                               sub,         &  ! subdomain
                               nn_sub,      &  ! number of nodes in the sub. 
                               nc_sub,      &  ! number of cells in the sub. 
                               nf_sub,      &  ! number of faces in the sub.
                               nbc_sub,     &  ! number of bnd. cells in sub
                               nbf_sub)        ! number of buffer cells in sub.
!------------------------------------------------------------------------------!
!   Writes file with connectivity data: name.cns                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, nn_sub, nc_sub, nf_sub,  &
                     nbc_sub,  nbf_sub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: b, c, s, n, lev, item, fu
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .cns file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, processor=sub, extension='.cns')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  write(fu) nn_sub
  write(fu) nc_sub + nbf_sub   ! new way: add buffer cells to cells
  write(fu) nbc_sub            ! number of boundary cells
  write(fu) nf_sub + nbf_sub
  write(fu) nbf_sub            ! number of buffer faces/cells
  write(fu) grid % n_bnd_cond  ! number of bounary conditions
  write(fu) grid % n_levels    ! number of multigrid levels

  !-------------------!
  !   Material name   !
  !-------------------!
  write(fu) grid % material % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, grid % n_bnd_cond
    write(fu) grid % bnd_cond % name(n)
  end do

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % cells_n_nodes(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % cells_n_nodes(buf_recv_ind(s))
  end do

  ! Cells' nodes
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      do n = 1, grid % cells_n_nodes(c)
        write(fu) grid % new_n(grid % cells_n(n,c))
      end do
    end if
  end do
  do s = 1, nbf_sub
    do n = 1, grid % cells_n_nodes(buf_recv_ind(s))
      write(fu) grid % new_n(grid % cells_n(n,buf_recv_ind(s)))
    end do
  end do

  ! Cells' processor ids
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % comm % proces(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % comm % proces(buf_recv_ind(s))
  end do
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % comm % proces(c)
    end if
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  do s = 1, grid % n_faces
    if(grid % new_f(s) .ne. 0) then
      write(fu) grid % faces_n_nodes(s)
    end if
  end do

  ! Faces' nodes
  do s = 1, grid % n_faces
    if(grid % new_f(s) .ne. 0) then
      do n = 1, grid % faces_n_nodes(s)
        write(fu) grid % new_n(grid % faces_n(n,s))
      end do
    end if
  end do

  ! Faces' cells
  do s = 1, grid % n_faces  ! OK, later chooses just faces with grid % new_f
    if( grid % new_f(s) > 0  .and.  grid % new_f(s) <= nf_sub ) then
      write(fu) grid % new_c(grid % faces_c(1,s)),  &
                grid % new_c(grid % faces_c(2,s))
    end if
  end do

  ! nbf_sub buffer faces (copy faces here, avoid them with buf_pos)
  do s = 1, nbf_sub
    if(buf_pos(s) > nc_sub) then     ! normal buffer (non-copy)
      write(fu) buf_send_ind(s),  &  ! new cell number
               buf_pos(s)            ! position in the buffer
    end if
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % bnd_cond % color(c)
    end if
  end do

  ! Boundary copy cells (this is a complete mess)
  do c = -1, -grid % n_bnd_cells, -1  ! OK, later chooses just ...
    if(grid % new_c(c) .ne. 0) then   ! ... cells with grid % new_c
      item = grid % new_c(grid % bnd_cond % copy_c(c))
      if(grid % bnd_cond % copy_c(c) .ne. 0) then
        if(grid % comm % proces(grid % bnd_cond % copy_c(c)) .ne. sub) then
          do b=1,nbf_sub
            if(buf_recv_ind(b) .eq. grid % bnd_cond % copy_c(c)) then
              print *, buf_pos(b) 
              print *, grid % xc(grid % bnd_cond % copy_c(c)),  &
                       grid % yc(grid % bnd_cond % copy_c(c)),  &
                       grid % zc(grid % bnd_cond % copy_c(c))
              item = -buf_pos(b) ! - sign, copy buffer
            end if
          end do
        endif
      endif
      write(fu) item
    end if
  end do

  !----------!
  !   Copy   !
  !----------!
  write(fu) grid % n_copy
  write(fu) ((grid % bnd_cond % copy_s(c,s), c = 1, 2), s = 1, grid % n_copy)

  !----------------------!
  !   Multigrid levels   !
  !----------------------!
  do lev = 1, grid % n_levels
    write(fu) grid % level(lev) % n_cells
    write(fu) grid % level(lev) % n_faces
  end do
  do lev = 1, grid % n_levels
    write(fu) (grid % level(lev) % cell(c),     c=1,grid % n_cells)
    write(fu) (grid % level(lev) % face(s),     s=1,grid % n_faces)
    write(fu) (grid % level(lev) % coarser_c(c),c=1,grid % level(lev) % n_cells)
    write(fu) (grid % level(lev) % faces_c(1,s),s=1,grid % level(lev) % n_faces)
    write(fu) (grid % level(lev) % faces_c(2,s),s=1,grid % level(lev) % n_faces)
  end do

  close(fu)

  end subroutine
