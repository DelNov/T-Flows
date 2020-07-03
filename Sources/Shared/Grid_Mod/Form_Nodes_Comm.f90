!==============================================================================!
  subroutine Grid_Mod_Form_Nodes_Comm(grid)
!------------------------------------------------------------------------------!
!   Find communication patterns for nodes from Process                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, c, c2, ln, n, s  ! counters
  integer :: max_n_procs         ! max number of processors surrounding a node

  integer :: nu, sub, p, cnt
  integer, allocatable :: nodes_cons(:)
  logical, allocatable :: active(:)     ! active nodes
  real,    allocatable :: xn_buff(:), yn_buff(:), zn_buff(:)
!==============================================================================!

  !------------------------------------------------------------------------!
  !   Allocate memory for number of processors in which the node resides   !
  !------------------------------------------------------------------------!
  n = grid % n_nodes
  allocate(grid % comm % nodes_n_procs(1:n))
  grid % comm % nodes_n_procs(:) = 0

  !-------------------------------------------!
  !   Allocate memory to store active nodes   !
  !-------------------------------------------!
  allocate(active(grid % n_nodes)); active(:) = .false.

  !-----------------------------------------------------!
  !   Find maximum number of processors for each node   !
  !-----------------------------------------------------!

  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      ! Increase possible number of processors ...
      ! ... surrounding the this node by one
      grid % comm % nodes_n_procs(n) =  &
      grid % comm % nodes_n_procs(n) + 1

      ! If current cell is not a buffer cell, make the node active
      if(c <= grid % n_cells - grid % comm % n_buff_cells) then
        active(n) = .true.
      end if
    end do
  end do

  max_n_procs = maxval(grid % comm % nodes_n_procs)

  ! Allocate memory for processors surrounding each node
  n = grid % n_nodes
  allocate(grid % comm % nodes_p(1:max_n_procs+1, 1:n))
  grid % comm % nodes_p(1:max_n_procs+1, 1:n) = HUGE_INT

  !---------------------------------------------------------------!
  !   Now you can really store the processors surrounding nodes   !
  !---------------------------------------------------------------!
  grid % comm % nodes_n_procs(:) = 0  ! re-initialize the cell count

  ! Inside cells
  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      if(active(n)) then

        ! Increase number of cells surrounding the this node by one ...
        grid % comm % nodes_n_procs(n) =  &
        grid % comm % nodes_n_procs(n) + 1

        ! ... and store the current cell and its processor
        grid % comm % nodes_p(grid % comm % nodes_n_procs(n), n)  &
                              = grid % comm % cell_proc(c)

      end if
    end do
  end do

  ! Find number of processors in each node by unique sort
  do n = 1, grid % n_nodes
    if(active(n)) then
      call Sort_Mod_Unique_Int(grid % comm % nodes_p(1:max_n_procs+1, n), nu)
      nu = nu - 1  ! remove the last HUGE_INT
      grid % comm % nodes_n_procs(n) = nu
    end if
  end do

  allocate(nodes_cons(n_proc)); nodes_cons(:) = 0

  ! Count node connections with other subdomains
  do n = 1, grid % n_nodes
    if(active(n)) then
      do p = 1, grid % comm % nodes_n_procs(n)
        sub = grid % comm % nodes_p(p, n)
        if(sub .ne. this_proc) then
          nodes_cons(sub) = nodes_cons(sub) + 1
        end if
      end do
    end if
  end do

  allocate(grid % comm % nodes_repl(n_proc))

  allocate(xn_buff(maxval(nodes_cons)));  xn_buff(:) = 0.0
  allocate(yn_buff(maxval(nodes_cons)));  yn_buff(:) = 0.0
  allocate(zn_buff(maxval(nodes_cons)));  zn_buff(:) = 0.0

  do sub = 1, n_proc
    n = nodes_cons(sub)
    grid % comm % nodes_repl(sub) % n_items = n
    if(n > 0) then
      allocate(grid % comm % nodes_repl(sub) % map(n));
      allocate(grid % comm % nodes_repl(sub) % i_buff(n));
      allocate(grid % comm % nodes_repl(sub) % l_buff(n));
      allocate(grid % comm % nodes_repl(sub) % r_buff(n));
    end if
  end do

  do sub = 1, n_proc
    cnt = 0
    if(grid % comm % nodes_repl(sub) % n_items > 0) then
      do n = 1, grid % n_nodes
        if( any(grid % comm % nodes_p(:,n) .eq. sub) ) then
          cnt = cnt + 1
          xn_buff(cnt) = grid % xn(n)
          yn_buff(cnt) = grid % yn(n)
          zn_buff(cnt) = grid % zn(n)
          grid % comm % nodes_repl(sub) % map(cnt) = n;
        end if
      end do

      ! Sort buffer cells by their coordinates
      call Sort_Mod_3_Real_Carry_Int(  &
                      xn_buff(1:cnt),  &
                      yn_buff(1:cnt),  &
                      zn_buff(1:cnt),  &
                      grid % comm % nodes_repl(sub) % map(1:cnt))

    end if
  end do

! do sub = 1, n_proc
!   do ln = 1, grid % comm % nodes_repl(sub) % n_items
!     n = grid % comm % nodes_repl(sub) % map(ln)
!     write(600 + this_proc*10 + sub, '(A,I5,A,I5,3F9.3)')           &
!                                 ' n=', n,                          &
!                                 ' g=', grid % comm % node_glo(n),  &
!                                                     grid % xn(n),  &
!                                                     grid % yn(n),  &
!                                                     grid % zn(n)
!   end do
! end do

  end subroutine
