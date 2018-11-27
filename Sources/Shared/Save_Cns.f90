!==============================================================================!
  subroutine Save_Cns(grid,        &
                      sub,         &  ! subdomain
                      nn_sub,      &  ! number of nodes in the sub. 
                      nc_sub,      &  ! number of cells in the sub. 
                      nf_sub,      &  ! number of faces in the sub.
                      nbc_sub,     &  ! number of bnd. cells in sub
                      nbf_sub)        ! number of buffer cells in sub.
!------------------------------------------------------------------------------!
!   Writes: name.cns                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, nn_sub, nc_sub, nf_sub,  &
                     nbc_sub,  nbf_sub
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c, s, n, c1, c2, count, var, lev
  integer              :: lower_bound, upper_bound
  character(len=80)    :: name_out
  integer, allocatable :: iwork(:,:)
  real, allocatable    :: work(:)
!==============================================================================!
!   The files name.cns and name.geo should merge into one file in some         !
!   of the future releases.                                                    !
!                                                                              !
!   sub     - subdomain number                                                 !
!   nn_sub  - number of nodes in subdomain                                     !
!   nc_sub  - number of cells in subdomain                                     !
!   nf_sub  - number of faces in subdomain, but without faces on buffer        !
!   nbc_sub - number of physicall boundary cells in subdomain                  !
!   nbf_sub - number of buffer boundary faces in subdomain                     !
!------------------------------------------------------------------------------!

  lower_bound = min(-nbf_sub, -grid % n_bnd_cells)
  upper_bound = max(grid % n_cells*8, grid % n_faces*4)

  allocate(iwork(lower_bound:upper_bound, 0:2));  iwork = 0
  allocate(work(grid % n_faces));                  work = 0.

  !----------------------!
  !                      !
  !   Create .cns file   !
  !                      !
  !----------------------!
  call Name_File( sub, name_out, '.cns' )
  open(9, file=name_out,form='unformatted', access='stream')
  write(*, *) '# Creating the file: ', trim(name_out)

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  write(9) nn_sub
  write(9) nc_sub + nbf_sub   ! new way: add buffer cells to cells
  write(9) nbc_sub            ! number of boundary cells
  write(9) nf_sub + nbf_sub
  write(9) nbf_sub            ! number of buffer faces/cells
  write(9) grid % n_bnd_cond  ! number of bounary conditions
  write(9) grid % n_levels    ! number of multigrid levels

  !-------------------!
  !   Material name   !
  !-------------------!
  write(9) grid % material % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = 1, grid % n_bnd_cond
    write(9) grid % bnd_cond % name(n)
  end do

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      count = count + 1
      iwork(count,1) = grid % cells_n_nodes(c)
    end if
  end do
  do s = 1, nbf_sub
    count = count + 1
    iwork(count,1) = grid % cells_n_nodes(buf_recv_ind(s))
  end do
  write(9) (iwork(c,1), c = 1, count)

  ! Cells' nodes
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      do n = 1, grid % cells_n_nodes(c)
        count = count + 1
        iwork(count,1) = grid % new_n(grid % cells_n(n,c))
      end do
    end if
  end do
  do s = 1, nbf_sub
    do n = 1, grid % cells_n_nodes(buf_recv_ind(s))
      count = count + 1
      iwork(count,1) = grid % new_n(grid % cells_n(n,buf_recv_ind(s)))
    end do
  end do
  write(9) (iwork(c,1), c = 1, count)

  ! Cells' processor ids
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) .ne. 0) then
      count = count + 1
      iwork(count,1) = grid % comm % proces(c)
    end if
  end do
  do s = 1, nbf_sub
    count = count + 1
    iwork(count,1) = grid % comm % proces(buf_recv_ind(s))
  end do
  write(9) (iwork(c,1), c = 1, count)

  ! Materials on boundary cells
  count = 0
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % new_c(c) .ne. 0) then
      count = count + 1
      iwork(count,1) = grid % comm % proces(c)
    end if
  end do
  write(9) (iwork(c,1), c = 1, count)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  count = 0
  do s = 1, grid % n_faces
    if(grid % new_f(s) .ne. 0) then
      count = count + 1
      iwork(count,1) = grid % faces_n_nodes(s)
    end if
  end do
  write(9) (iwork(s,1), s = 1, count)

  ! Faces' nodes
  count = 0
  do s = 1, grid % n_faces
    if(grid % new_f(s) .ne. 0) then
      do n = 1, grid % faces_n_nodes(s)
        count = count + 1
        iwork(count,1) = grid % new_n(grid % faces_n(n,s))
      end do
    end if
  end do
  write(9) (iwork(s,1), s = 1, count)

  count = 0

  ! nf_sub physical faces
  do s = 1, grid % n_faces  ! OK, later chooses just faces with grid % new_f
    if( grid % new_f(s) > 0  .and.  grid % new_f(s) <= nf_sub ) then
      count = count + 1
      iwork(count,1) = grid % new_c(grid % faces_c(1,s))
      iwork(count,2) = grid % new_c(grid % faces_c(2,s))
    end if
  end do

  ! nbf_sub buffer faces (copy faces here, avoid them with buf_pos)
  do s = 1, nbf_sub
    if(buf_pos(s) > nc_sub) then        ! normal buffer (non-copy)
      count = count + 1
      iwork(count,1) = buf_send_ind(s)  ! new cell number
      iwork(count,2) = buf_pos(s)       ! position in the buffer
    end if
  end do

  write(9) (iwork(s,1), s = 1, count)
  write(9) (iwork(s,2), s = 1, count)

  !--------------!
  !   Boundary   !
  !--------------!
  count = 0 ! count goes to negative

  ! nbc_sub physical boundary cells
  do c = -1, -grid % n_bnd_cells, -1  ! OK, later chooses just cells with grid % new_c
    if(grid % new_c(c) .ne. 0) then
      count=count-1
      ! nekad bio i: grid % new_c(c)
      iwork(count,1) = grid % bnd_cond % color(c)
      iwork(count,2) = grid % new_c(grid % bnd_cond % copy_c(c))
      if(grid % bnd_cond % copy_c(c) .ne. 0) then
        if(grid % comm % proces(grid % bnd_cond % copy_c(c)) .ne. sub) then
          do b=1,nbf_sub
            if(buf_recv_ind(b) .eq. grid % bnd_cond % copy_c(c)) then
              print *, buf_pos(b) 
              print *, grid % xc(grid % bnd_cond % copy_c(c)),  &
                       grid % yc(grid % bnd_cond % copy_c(c)),  &
                       grid % zc(grid % bnd_cond % copy_c(c))
              iwork(count,2)=-buf_pos(b) ! - sign, copy buffer
            end if
          end do
        endif
      endif
    end if
  end do 

  write(9) (iwork(c,1), c = -1, count, -1)
  write(9) (iwork(c,2), c = -1, count, -1)

  !----------!
  !   Copy   !
  !----------!
  count = 0
  do s = 1, grid % n_copy
    count = count + 1
    iwork(count,1) = grid % bnd_cond % copy_s(1,s)
    iwork(count,2) = grid % bnd_cond % copy_s(2,s)
  end do

  write(9) count 
  write(9) (iwork(c,1), c = 1, count)
  write(9) (iwork(c,2), c = 1, count)

  !----------------------!
  !   Multigrid levels   !
  !----------------------!
  do lev = 1, grid % n_levels
    write(9) grid % level(lev) % n_cells
    write(9) grid % level(lev) % n_faces
  end do
  do lev = 1, grid % n_levels
    write(9) (grid % level(lev) % cell(c),      c=1,grid % n_cells)
    write(9) (grid % level(lev) % face(s),      s=1,grid % n_faces)
    write(9) (grid % level(lev) % faces_c(1,s), s=1,grid % level(lev) % n_faces)
    write(9) (grid % level(lev) % faces_c(2,s), s=1,grid % level(lev) % n_faces)
  end do
  close(9)

  deallocate (iwork)
  deallocate (work)

  end subroutine
