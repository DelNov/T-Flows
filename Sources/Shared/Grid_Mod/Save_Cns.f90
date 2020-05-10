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
  integer           :: b, c, s, n, lev, item, fu, mu
  character(len=80) :: name_out
  character(len=80) :: name_map
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .map file   !
  !                      !
  !----------------------!
  if(sub > 0) then
    call File_Mod_Set_Name(name_map, processor=sub, extension='.map')
    call File_Mod_Open_File_For_Writing(name_map, mu)
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        write(mu, '(i9)') grid % comm % cell_glo(c)
      end if
    end do
    do c = -grid % n_bnd_cells, -1
      if(grid % comm % cell_proc(c) .eq. sub) then
        write(mu, '(i9)') grid % comm % cell_glo(c)
      end if
    end do
    close(mu)
  end if   ! through subdomains

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
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % cells_n_nodes(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % cells_n_nodes(buf_recv_ind(s))
  end do

  ! Cells' nodes
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .eq. sub) then
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
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % comm % cell_proc(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % comm % cell_proc(buf_recv_ind(s))
  end do
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % comm % cell_proc(c)
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

  ! nbf_sub buffer faces and the new numbers of cells surrounding them
  do s = 1, nbf_sub
    write(fu) buf_send_ind(s),  &
              nc_sub + s
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % comm % cell_proc(c) .eq. sub) then
      write(fu) grid % bnd_cond % color(c)
    end if
  end do

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
