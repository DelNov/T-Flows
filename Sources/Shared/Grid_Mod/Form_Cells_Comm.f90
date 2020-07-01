!==============================================================================!
  subroutine Grid_Mod_Form_Cells_Comm(grid)
!------------------------------------------------------------------------------!
!   Find communication patterns for cells from Process                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: sub, ln, ms, mr, n, nn, sh
  character(len=80)    :: name_in
  integer              :: c, c1, c2, s, buf_cnt
  real,    allocatable :: xf_buff(:), yf_buff(:), zf_buff(:)
!==============================================================================!

  if(n_proc < 2) return

  !--------------------------------------------------!
  !   Allocate memory for send and receive buffers   !
  !--------------------------------------------------!
  allocate(grid % comm % cells_send(n_proc))
  allocate(grid % comm % cells_recv(n_proc))

  !------------------------!
  !   Count buffer cells   !
  !------------------------!
  grid % comm % n_buff_cells = 0
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .ne. this_proc) then
      grid % comm % n_buff_cells =    &
      grid % comm % n_buff_cells + 1
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
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 > 0) then
          if(grid % comm % cell_proc(c1) .eq. this_proc .and.  &
             grid % comm % cell_proc(c2) .eq. sub) then
            do ln = 1, grid % faces_n_nodes(s)
              grid % new_n(grid % faces_n(ln, s)) = -1
            end do
          end if
        end if
      end do

      ! Spread info to twin nodes
      do s = 1, grid % n_faces
        if(grid % faces_s(s) > 0) then
          c1 = grid % faces_c(1, s)
          c2 = grid % faces_c(2, s)
          if(grid % comm % cell_proc(c1) .eq. this_proc .and.  &
             grid % comm % cell_proc(c2) .eq. sub       .or.   &
             grid % comm % cell_proc(c2) .eq. this_proc .and.  &
             grid % comm % cell_proc(c1) .eq. sub) then
            nn = grid % faces_n_nodes(s)
            sh = grid % faces_s(s)
            grid % new_n(grid % faces_n(1:nn, s )) = -1
            grid % new_n(grid % faces_n(1:nn, sh)) = -1
          end if
        end if
      end do

      ms = 0
      mr = 0
      do c = 1, grid % n_cells
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
      do c = 1, grid % n_cells
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

  end subroutine
