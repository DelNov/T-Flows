!==============================================================================!
  subroutine Grid_Mod_Form_Cells_Comm(grid)
!------------------------------------------------------------------------------!
!   Find communication patterns for cells from Process                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: sub, i_nod, j_nod, ms, mr, n, ni, nj, sh
  character(SL) :: name_in
  integer       :: c, c1, c2, s, n_buff_faces
  real          :: xi, yi, zi, xj, yj, zj
!==============================================================================!

  if(n_proc < 2) return

  !--------------------------------------------------!
  !   Allocate memory for send and receive buffers   !
  !--------------------------------------------------!
  allocate(grid % comm % cells_send(n_proc))
  allocate(grid % comm % cells_recv(n_proc))

  !-------------------------------!
  !   Count buffer cells inside   !
  !-------------------------------!
  grid % comm % n_buff_cells = 0
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .ne. this_proc) then
      grid % comm % n_buff_cells =    &
      grid % comm % n_buff_cells + 1
    end if
  end do

  !------------------------!
  !   Count buffer faces   !
  !------------------------!
  n_buff_faces = 0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(grid % comm % cell_proc(c1) .ne. this_proc .and.  &
       grid % comm % cell_proc(c2) .ne. this_proc) then
      n_buff_faces = n_buff_faces + 1
    end if
  end do

  allocate(grid % new_n(grid % n_nodes))

  !--------------------------------------------------!
  !   Calculate size of buffers between subdomains   !
  !--------------------------------------------------!
  do sub = 1, n_proc

    ! Initialize buffer size to zero
    grid % comm % cells_send(sub) % n_items = 0
    grid % comm % cells_recv(sub) % n_items = 0

    if(sub .ne. this_proc) then

      ! Mark the nodes at subdomain interface
      grid % new_n(:) = 0
      do c = 1, grid % n_cells
        if(grid % comm % cell_proc(c) .eq. sub) then
          do i_nod = 1, grid % cells_n_nodes(c)
            grid % new_n(grid % cells_n(i_nod, c)) = -1
          end do
        end if
      end do

      ! Spread info to twin nodes
      ! (This is not super elegant, but seems to work)
      do s = 1, grid % n_faces
        sh = grid % faces_s(s)
        if(sh > 0) then
          do i_nod = 1, grid % faces_n_nodes(s)
            ni = grid % faces_n(i_nod, s)
            if(grid % new_n(ni) .eq. -1) then
              xi = grid % xn(ni);  yi = grid % yn(ni); zi = grid % zn(ni)
              do j_nod = 1, grid % faces_n_nodes(sh)
                nj = grid % faces_n(j_nod, sh)
                xj = grid % xn(nj);  yj = grid % yn(nj);  zj = grid % zn(nj)
                if(grid % per_x > TINY) then
                  if(Math_Mod_Distance(0., yi, zi, 0., yj, zj) < NANO) then
                    grid % new_n(nj) = -1
                  end if
                end if
                if(grid % per_y > TINY) then
                  if(Math_Mod_Distance(xi, 0., zi, xj, 0., zj) < NANO) then
                    grid % new_n(nj) = -1
                  end if
                end if
                if(grid % per_z > TINY) then
                  if(Math_Mod_Distance(xi, yi, 0., xj, yj, 0.) < NANO) then
                    grid % new_n(nj) = -1
                  end if
                end if
              end do
            end if
          end do
          do j_nod = 1, grid % faces_n_nodes(sh)
            nj = grid % faces_n(j_nod, sh)
            if(grid % new_n(nj) .eq. -1) then
              xj = grid % xn(nj);  yj = grid % yn(nj);  zj = grid % zn(nj)
              do i_nod = 1, grid % faces_n_nodes(s)
                ni = grid % faces_n(i_nod, s)
                xi = grid % xn(ni);  yi = grid % yn(ni); zi = grid % zn(ni)
                if(grid % per_x > TINY) then
                  if(Math_Mod_Distance(0., yi, zi, 0., yj, zj) < NANO) then
                    grid % new_n(ni) = -1
                  end if
                end if
                if(grid % per_y > TINY) then
                  if(Math_Mod_Distance(xi, 0., zi, xj, 0., zj) < NANO) then
                    grid % new_n(ni) = -1
                  end if
                end if
                if(grid % per_z > TINY) then
                  if(Math_Mod_Distance(xi, yi, 0., xj, yj, 0.) < NANO) then
                    grid % new_n(ni) = -1
                  end if
                end if
              end do
            end if
          end do
        end if  ! sh > 0
      end do

      ms = 0
      mr = 0
      do c = -grid % n_bnd_cells, grid % n_cells
        if(grid % comm % cell_proc(c) .eq. this_proc) then
          n = grid % cells_n_nodes(c)
          if( any( grid % new_n(grid % cells_n(1:n,c)) .eq. -1) ) then
            ms = ms + 1
          end if
        end if
        if(grid % comm % cell_proc(c) .eq. sub) then
          n = grid % cells_n_nodes(c)
          if( any( grid % new_n(grid % cells_n(1:n,c)) .eq. -1) ) then
            mr = mr + 1
          end if
        end if
      end do

      if(ms > 0) then
        allocate(grid % comm % cells_send(sub) % map(ms));
        allocate(grid % comm % cells_send(sub) % i_buff(ms));
        allocate(grid % comm % cells_send(sub) % l_buff(ms));
        allocate(grid % comm % cells_send(sub) % r_buff(ms));
      end if
      if(mr > 0) then
        allocate(grid % comm % cells_recv(sub) % map(mr));
        allocate(grid % comm % cells_recv(sub) % i_buff(mr));
        allocate(grid % comm % cells_recv(sub) % l_buff(mr));
        allocate(grid % comm % cells_recv(sub) % r_buff(mr));
      end if

      ms = 0
      mr = 0
      do c = -grid % n_bnd_cells, grid % n_cells
        if(grid % comm % cell_proc(c) .eq. this_proc) then
          n = grid % cells_n_nodes(c)
          if( any( grid % new_n(grid % cells_n(1:n,c)) .eq. -1) ) then
            ms = ms + 1
            grid % comm % cells_send(sub) % map(ms) = c
          end if
        end if
        if(grid % comm % cell_proc(c) .eq. sub) then
          n = grid % cells_n_nodes(c)
          if( any( grid % new_n(grid % cells_n(1:n,c)) .eq. -1) ) then
            mr = mr + 1
            grid % comm % cells_recv(sub) % map(mr) = c
          end if
        end if
      end do

      ! Store final buffer lengths
      grid % comm % cells_send(sub) % n_items = ms
      grid % comm % cells_recv(sub) % n_items = mr

    end if

  end do

  grid % n_faces = grid % n_faces - n_buff_faces

  end subroutine
