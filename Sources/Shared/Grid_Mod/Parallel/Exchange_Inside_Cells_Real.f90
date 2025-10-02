!==============================================================================!
  subroutine Exchange_Inside_Cells_Real(Grid, phi, caller)
!------------------------------------------------------------------------------!
!>  The subroutine's purpose is to synchronize the values of an integer array
!>  (phi) that is associated with the cells of the computational grid (Grid).
!>  This array is distributed across various processors, and the subroutine
!>  ensures that each processor receives the necessary data from its
!>  neighboring processors.  In 2021, a neat extension was introduced to this
!>  procedure.  One can send an optional parameter (caller), a string, with
!>  name and line of the caller procedure.  (Such as: 'Compute_Momentum_347'.)
!>  If that parameter is sent, this procedure will check if the call was
!>  needed or redunant and the info is stored in file 'exchange_cells_real.log'.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparing send buffers:                                                  !
!     - Iterates over cells that need to send data. It maps values from phi    !
!       to a send buffer (r_buff).                                             !
!   * Data exchange:                                                           !
!     - Performs a send-receive operation across processors. Sends data from   !
!       r_buff and simultaneously receives data into r_buff of cells_recv.     !
!   * Updating local data:                                                     !
!     - Updates phi array with new values received, ensuring consistency       !
!       across processors.                                                     !
!   * Optional efficiency check:                                               !
!     - If the caller parameter is provided, the subroutine logs whether       !
!       the call was necessary. This is done by comparing old and new values   !
!       in the receive buffer, and logging the result in                       !
!       'exchange_cells_real.log'.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)           :: Grid    !! computational grid
  real, intent(inout)        :: phi(1:Grid % n_cells)
    !! real array defined over boundary and internal cells to be exchanged
  character(len=*), optional :: caller  !! string to log call's necessity
!-----------------------------------[Locals]-----------------------------------!
  integer       :: ln, c1, c2, len_s, len_r, sub
  integer       :: needed
  integer, save :: fu
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES - 1
!==============================================================================!

  !-------------------------------------------------------------------!
  !   Optional parameter was sent; store old values for later check   !
  !-------------------------------------------------------------------!
  if(present(caller)) then

    if(First_Proc()) then
      call File % Append_For_Writing_Ascii(  &
                  'exchange_cells_real.log', fu, This_Proc())
    end if

    do sub = 1, N_Procs()
      len_r = Grid % Comm % cells_recv(sub) % n_items
      do ln = 1, len_r
        c2 = Grid % Comm % cells_recv(sub) % map(ln)
        Grid % Comm % cells_recv(sub) % o_buff(ln) = phi(c2)
      end do
    end do
  end if

  !----------------------------------!
  !   Usual course of this routine   !
  !----------------------------------!

  ! Fill the buffers with new values
  do sub = 1, N_Procs()
    len_s = Grid % Comm % cells_send(sub) % n_items
    do ln = 1, len_s
      c1 = Grid % Comm % cells_send(sub) % map(ln)
      Grid % Comm % cells_send(sub) % r_buff(ln) = phi(c1)
    end do
  end do

  ! Exchange the values
  do sub = 1, N_Procs()
    len_s = Grid % Comm % cells_send(sub) % n_items
    len_r = Grid % Comm % cells_recv(sub) % n_items
    if(len_s + len_r > 0) then
      call Grid % Comm % Sendrecv_Real_Arrays(      &
        len_s,                                      &  ! sending length
        Grid % Comm % cells_send(sub) % r_buff(1),  &  ! array to be sent
        len_r,                                      &  ! receiving length
        Grid % Comm % cells_recv(sub) % r_buff(1),  &  ! array to be received
        sub)                                           ! destination processor
    end if
  end do

  ! Fill the buffers with new values
  do sub = 1, N_Procs()
    len_r = Grid % Comm % cells_recv(sub) % n_items
    do ln = 1, len_r
      c2 = Grid % Comm % cells_recv(sub) % map(ln)
      phi(c2) = Grid % Comm % cells_recv(sub) % r_buff(ln)
    end do
  end do

  !-------------------------------------------------!
  !   Optional parameter was sent; check new with   !
  !    old values to see if this call was needed    !
  !-------------------------------------------------!
  if(present(caller)) then

    needed = NO
    do sub = 1, N_Procs()
      len_r = Grid % Comm % cells_recv(sub) % n_items
      do ln = 1, len_r
        if( Grid % Comm % cells_recv(sub) % r_buff(ln) .ne. &
            Grid % Comm % cells_recv(sub) % o_buff(ln) ) then
          needed = YES
        end if
      end do
    end do

    ! What if it is needed on some other processor
    call Global % Max_Int(needed)

    if(First_Proc()) then
      if(needed .eq. YES) then
        write(fu,'(a)') '# Call from: ' // caller // ' was needed!'
      else
        write(fu,'(a)') '# Call from: ' // caller // ' was redundant!'
      end if
      close(fu)
    end if

  end if

  end subroutine
