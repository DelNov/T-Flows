!==============================================================================!
  subroutine Form_Cells_Comm(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine plays a critical role in T-Flows for parallel processing.
!>  It ensures that each sub-domain correctly identifies and communicates with
!>  buffer cells from other sub-domains, maintaining the integrity and of the
!>  simulation process in a parallel computing environment.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Buffer cell analysis: Identifies buffer cells in each sub-domain and     !
!     records their global numbers and associated sub-domain numbers.          !
!   * Global information sharing: Distributes the buffer cell information      !
!     among all processors via a global sum over integer arrays, forming the   !
!     basis for buffer creation.                                               !
!   * Buffer formation: Based on the gathered information, it forms send and   !
!     receive buffers for each processor, identifying which cells need to send !
!     or receive data.                                                         !
!   * Sorting for robustness: Implements sorting mechanisms (send_sort() and   !
!     recv_sort()) to ensure correct ordering in send and receive buffers.     !
!     This is vital for consistency in data exchange across processors.        !
!   * Validation and debugging: Includes debugging provisions for verifying    !
!     the correct formation of communication patterns and ensuring that the    !
!     setup for parallel processing is accurate.                               !
!   * Error handling: In case of discrepancies (e.g., mismatched global        !
!     numbers in buffer cells), the subroutine triggers a critical error,      !
!     indicating potential setup issues in the parallel run configuration.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! computational grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer              :: sub, ms, mr, cnt, n_fail, i, j, i_cel, j_cel, i_nod
  integer              :: c, c1, c2, s, n_buff_faces, n_max_buff_cells, n
  integer              :: n_cells_near_buffers
  integer, allocatable :: cells_near_buffers(:)
  integer, allocatable :: send_cells(:), recv_cells(:)
  integer, allocatable :: send_buff_cnt(:,:), recv_buff_cnt(:,:)
  integer, allocatable :: need_cell(:,:), from_proc(:,:)
  integer, allocatable :: send_sort(:), recv_sort(:)
  integer, allocatable :: glo(:)  ! used for checking the communication
  logical, allocatable :: node_in_buffer(:)
!==============================================================================!

  if(Sequential_Run()) return

  if(First_Proc()) print '(a)', ' # Forming cells communication patterns ...'

  call Profiler % Start('Form_Cells_Com')

  ! Allocate memory for locally used arrays
  allocate(send_buff_cnt(N_Procs(), N_Procs()))
  allocate(recv_buff_cnt(N_Procs(), N_Procs()))
  allocate(send_cells(-Grid % n_bnd_cells:Grid % n_cells))
  allocate(recv_cells(-Grid % n_bnd_cells:Grid % n_cells))

  ! Find which cells are near buffers, it speeds up the procedure quite a bit
  allocate(node_in_buffer(Grid % n_nodes));  node_in_buffer(:) = .false.
  do c = Cells_In_Buffers()
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      node_in_buffer(n) = .true.
    end do
  end do
  n_cells_near_buffers = 0
  do c = Cells_In_Domain()
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      if(node_in_buffer(n)) then
        n_cells_near_buffers = n_cells_near_buffers + 1
        exit
      end if
    end do
  end do
  allocate(cells_near_buffers(n_cells_near_buffers))
  n_cells_near_buffers = 0
  do c = Cells_In_Domain()
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      if(node_in_buffer(n)) then
        n_cells_near_buffers = n_cells_near_buffers + 1
        cells_near_buffers(n_cells_near_buffers) = c
        exit
      end if
    end do
  end do

  !-------------------------------!
  !   Count buffer cells inside   !
  !-------------------------------!
  Grid % Comm % n_buff_cells = 0
  do c = Cells_In_Domain()
    Assert(Cell_In_This_Proc(c))
  end do
  do c = Cells_In_Buffers()
    Assert(.not. Cell_In_This_Proc(c))
    Grid % Comm % n_buff_cells = Grid % Comm % n_buff_cells + 1
  end do

  !------------------------!
  !   Count buffer faces   !
  !------------------------!
  n_buff_faces = 0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c1 .eq. 0 .and. c2 .eq. 0) then
      n_buff_faces = n_buff_faces + 1
    end if
  end do

  !--------------------------------------------------------!
  !   Form the recv_buff_cnt matrix and recv_cells array   !
  !--------------------------------------------------------!
  recv_buff_cnt(:,:)   = 0
  recv_cells(:) = 0
  do c = Cells_In_Buffers()
    sub = Grid % Comm % cell_proc(c)
    recv_buff_cnt(This_Proc(), sub) = recv_buff_cnt(This_Proc(), sub) + 1
    recv_cells(c) = sub
  end do

  if(DEBUG) then  ! a)
    do sub=1, N_Procs()
      if(sub.ne.This_Proc() .and. recv_buff_cnt(This_Proc(), sub).gt.0) then
        write(100*This_Proc()+sub, *)                                          &
                                    '#====================================' // &
                                    '====================================='
        write(100*This_Proc()+sub, '(a,i7,a,i7)')                    &
              ' # There are        ', Grid % Comm % n_buff_cells,  &
              ' buffer cells in processor ', This_Proc()
        write(100*This_Proc()+sub, *)                                          &
                                    '#------------------------------------' // &
                                    '-------------------------------------'
        write(100*This_Proc()+sub, '(a,i7,a,i7)')                       &
              ' a)  It needs       ', recv_buff_cnt(This_Proc(), sub),  &
              ' cells from processors     ', sub
      end if
    end do
  end if

  !--------------------------------------------------!
  !   Allocate memory for matrix with needed cells   !
  !--------------------------------------------------!
  n_max_buff_cells = Grid % Comm % n_buff_cells
  call Global % Max_Int(n_max_buff_cells)
  allocate(need_cell(n_max_buff_cells, N_Procs()));  need_cell(:,:) = 0
  allocate(from_proc(n_max_buff_cells, N_Procs()));  from_proc(:,:) = 0

  !-------------------------------------------!
  !   Store global number of cells you need   !
  !    and also from which processor it is    !
  !-------------------------------------------!
  i_cel = 0
  do c = Cells_In_Buffers()
    Assert(.not. Cell_In_This_Proc(c))
    i_cel = i_cel + 1
    need_cell(i_cel, This_Proc()) = Grid % Comm % cell_glo(c)
    from_proc(i_cel, This_Proc()) = Grid % Comm % cell_proc(c)
  end do

  !----------------------------------------------!
  !   Inform all processors about needed cells   !
  !   and from which processor are they needed   !
  !----------------------------------------------!
  need_cell = reshape(need_cell, (/n_max_buff_cells*N_Procs(), 1/))
  from_proc = reshape(from_proc, (/n_max_buff_cells*N_Procs(), 1/))
  call Global % Sum_Int_Array(n_max_buff_cells*N_Procs(), need_cell)
  call Global % Sum_Int_Array(n_max_buff_cells*N_Procs(), from_proc)
  need_cell = reshape(need_cell, (/n_max_buff_cells, N_Procs()/))
  from_proc = reshape(from_proc, (/n_max_buff_cells, N_Procs()/))

  !---------------------------------------------------------------------!
  !   Form send_buff_cnt from the information in the matrix from_proc   !
  !---------------------------------------------------------------------!
  send_buff_cnt(:,:) = 0
  do sub = 1, N_Procs()
    if(sub .ne. This_Proc()) then
      do i_cel = 1, n_max_buff_cells
        if(from_proc(i_cel, sub) .eq. This_Proc()) then
          send_buff_cnt(This_Proc(), sub) = send_buff_cnt(This_Proc(), sub) + 1
        end if
      end do
      if(DEBUG) then  ! b)
        if(send_buff_cnt(This_Proc(), sub).gt.0) then
          write(100*This_Proc() + sub, '(2(a,i7))')                           &
                    ' b)  It should send ', send_buff_cnt(This_Proc(), sub),  &
                    ' cells to processor        ', sub
        end if
      end if
    end if
  end do

  Assert(send_buff_cnt(This_Proc(), This_Proc()) .eq. 0)

  if(DEBUG) then
    send_buff_cnt = reshape(send_buff_cnt, (/N_Procs()*N_Procs(), 1/))
    recv_buff_cnt = reshape(recv_buff_cnt, (/N_Procs()*N_Procs(), 1/))
    call Global % Sum_Int_Array(N_Procs()*N_Procs(), send_buff_cnt)
    call Global % Sum_Int_Array(N_Procs()*N_Procs(), recv_buff_cnt)
    send_buff_cnt = reshape(send_buff_cnt, (/N_Procs(),N_Procs()/))
    recv_buff_cnt = reshape(recv_buff_cnt, (/N_Procs(),N_Procs()/))

    do i = 1, N_Procs()
      do j = 1, N_Procs()
        Assert(send_buff_cnt(i,j) .eq. recv_buff_cnt(j,i))
      end do
    end do
  end if

  !--------------------------------------------------!
  !   Allocate memory for send and receive buffers   !
  !--------------------------------------------------!
  allocate(Grid % Comm % cells_send(N_Procs()))
  allocate(Grid % Comm % cells_recv(N_Procs()))

  call Profiler % Start('Form_Cells_Com, Form Send/Recv Buffers')

  !-----------------------------------!
  !                                   !
  !   Form send and receive buffers   !
  !                                   !
  !-----------------------------------!
  do sub = 1, N_Procs()

    ! Initialize buffer size to zero
    Grid % Comm % cells_send(sub) % n_items = 0
    Grid % Comm % cells_recv(sub) % n_items = 0

    if(sub .ne. This_Proc()) then

      !---------------------------------!
      !   Allocate memory for buffers   !
      !---------------------------------!
      ms = send_buff_cnt(This_Proc(), sub)
      mr = recv_buff_cnt(This_Proc(), sub)

      if(ms > 0) then
        allocate(Grid % Comm % cells_send(sub) % map(ms));
        allocate(Grid % Comm % cells_send(sub) % i_buff(ms));
        allocate(Grid % Comm % cells_send(sub) % l_buff(ms));
        allocate(Grid % Comm % cells_send(sub) % r_buff(ms));
        allocate(send_sort(ms));
      end if

      if(mr > 0) then
        allocate(Grid % Comm % cells_recv(sub) % map(mr));
        allocate(Grid % Comm % cells_recv(sub) % i_buff(mr));
        allocate(Grid % Comm % cells_recv(sub) % l_buff(mr));
        allocate(Grid % Comm % cells_recv(sub) % r_buff(mr));
        allocate(Grid % Comm % cells_recv(sub) % o_buff(mr));
        allocate(recv_sort(mr));
      end if

      !------------------------------------------!
      !   Form send_cells from the information   !
      !   in matrices from_proc and need_cell    !
      !------------------------------------------!

      cnt = 0
      send_cells(:) = 0
      do i_cel = 1, n_max_buff_cells
        if(from_proc(i_cel, sub) .eq. This_Proc()) then
          do j_cel = 1, n_cells_near_buffers
            c = cells_near_buffers(j_cel)
            if(Grid % Comm % cell_glo(c) .eq. need_cell(i_cel, sub)) then
              send_cells(c) = sub                                  ! identify
              cnt = cnt + 1
              goto 1
            end if
          end do
        end if
1       continue
      end do
      if(DEBUG) then  ! c)
        if(send_buff_cnt(This_Proc(), sub).gt.0) then
          write(100*This_Proc() + sub, '(2(a,i7))')    &
                       ' c)  It did find    ', cnt,    &
                       ' cells to send to processor', sub
        end if
      end if

      ! Browse through cells in domain to find what you have to send
      ms = 0
      do c = Cells_In_Domain()
        if(send_cells(c) .eq. sub) then
          ms = ms + 1
          Grid % Comm % cells_send(sub) % map(ms) = c
          send_sort(ms) = Grid % Comm % cell_glo(c)    ! store sorting criterion
        end if
      end do

      ! Important link for robustness: insures that the sent cells in this
      ! domain are ordered in the same way as receive cells in other domain
      if(ms > 0) then
        call Sort % Int_Carry_Int(send_sort,  &
                                  Grid % Comm % cells_send(sub) % map)
      end if

      ! Browse through cells in buffers to find out what you need to receive
      mr = 0
      do c = Cells_In_Buffers()
        if(recv_cells(c) .eq. sub) then
          mr = mr + 1
          Grid % Comm % cells_recv(sub) % map(mr) = c
          recv_sort(mr) = Grid % Comm % cell_glo(c)    ! store sorting criterion
        end if
      end do

      ! Important link for robustness: insures that the receive cells in this
      ! domain are ordered in the same way as sent cells in other domain.
      ! (Yet, this is not engaged yet since the cells in buffers has not been
      !  fiddled with like the cells in domain due to OpenMP threading.)
      ! if(mr > 0) then
      !   call Sort % Int_Carry_Int(recv_sort,  &
      !                             Grid % Comm % cells_recv(sub) % map)
      ! end if

      if(DEBUG) then  !  d)
        if(recv_buff_cnt(This_Proc(), sub).gt.0 .and.  &
           send_buff_cnt(This_Proc(), sub).gt.0) then
          write(100*This_Proc()+sub, '(a,i0.0,a,i0.0,a,i0.0,a,i0.0)')  &
              ' d)  send/recv (', This_Proc(), '/', sub, ') =  ', ms, ' / ', mr
        end if
      end if

      ! Store final buffer lengths
      Grid % Comm % cells_send(sub) % n_items = ms
      Grid % Comm % cells_recv(sub) % n_items = mr

      ! Free sorting arrays for the next iteration
      if(ms > 0) deallocate(send_sort)
      if(mr > 0) deallocate(recv_sort)

    end if  ! sub .ne. This_Proc()
  end do    ! sub

  call Profiler % Stop('Form_Cells_Com, Form Send/Recv Buffers')

  !-------------------------------------!
  !   Avoid faces in the buffers only   !
  !   (This also leaves shadows out)    !
  !-------------------------------------!
  Grid % n_faces = Grid % n_faces - n_buff_faces

  !------------------------------------------!
  !   Check if communication patterns work   !
  !------------------------------------------!
  allocate(glo(-Grid % n_bnd_cells : Grid % n_cells));  glo(:) = 0

  do c = Cells_In_Domain()
    glo(c) = Grid % Comm % cell_glo(c)
  end do

  call Grid % Exchange_Cells_Int(glo)

  n_fail = 0
  do c = Cells_In_Buffers()
    if(.not. glo(c) .eq. Grid % Comm % cell_glo(c)) then
      n_fail = n_fail + 1
    end if
  end do
  call Global % Sum_Int(n_fail)
  if(n_fail .gt. 0) then
    call Message % Error(55,                                                  &
                         'Ouch, this hurts. Formation of communication '  //  &
                         'patterns has failed. Hopefully, you just ran '  //  &
                         'the simulation on fewer processor than there '  //  &
                         'are subdomains. \n \n This error is critical.',     &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  if(First_Proc()) print '(a)', ' # Done !'

  call Profiler % Stop('Form_Cells_Com')

  end subroutine
