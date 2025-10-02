!==============================================================================!
  subroutine Exchange_Cells_Int(Grid, phi)
!------------------------------------------------------------------------------!
!>  The subroutine's purpose is to synchronize the values of an integer array
!>  (phi) that is associated with the cells of the computational grid (Grid).
!>  This array is distributed across various processors, and the subroutine
!>  ensures that each processor receives the necessary data from its
!>  neighboring processors.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparing send buffers:                                                  !
!     - For each processor (sub), the subroutine iterates over the cells that  !
!       need to send data (cells_send). It maps the relevant values from the   !
!       phi array to a send buffer (i_buff).                                   !
!   * Data exchange:                                                           !
!     - The subroutine performs a send-receive operation for each processor.   !
!       The data from the send buffer (i_buff of cells_send) is sent to the    !
!       corresponding processor, and simultaneously, data is received from the !
!       other processors into the receive buffer (i_buff of cells_recv).       !
!   * Updating local data:                                                     !
!     - After receiving data, the subroutine updates the phi array with the    !
!       new values from the receive buffer. This step ensures that each        !
!       processor has the updated data from its neighboring processors.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)        :: Grid  !! computational grid
  integer, intent(inout)  :: phi(-Grid % n_bnd_cells:Grid % n_cells)
    !! integer array defined over bounardy and internal cells to be exchanged
!-----------------------------------[Locals]-----------------------------------!
  integer :: ln, c1, c2, len_s, len_r, sub
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, N_Procs()
    len_s = Grid % Comm % cells_send(sub) % n_items
    do ln = 1, len_s
      c1 = Grid % Comm % cells_send(sub) % map(ln)
      Grid % Comm % cells_send(sub) % i_buff(ln) = phi(c1)
    end do
  end do

  ! Exchange the values
  do sub = 1, N_Procs()
    len_s = Grid % Comm % cells_send(sub) % n_items
    len_r = Grid % Comm % cells_recv(sub) % n_items
    if(len_s + len_r > 0) then
      call Grid % Comm % Sendrecv_Int_Arrays(       &
        len_s,                                      &  ! sending length
        Grid % Comm % cells_send(sub) % i_buff(1),  &  ! array to be sent
        len_r,                                      &  ! receiving length
        Grid % Comm % cells_recv(sub) % i_buff(1),  &  ! array to be received
        sub)                                           ! destination processor
    end if
  end do

  ! Fill the buffers with new values
  do sub = 1, N_Procs()
    len_r = Grid % Comm % cells_recv(sub) % n_items
    do ln = 1, len_r
      c2 = Grid % Comm % cells_recv(sub) % map(ln)
      phi(c2) = Grid % Comm % cells_recv(sub) % i_buff(ln)
    end do
  end do

  end subroutine
