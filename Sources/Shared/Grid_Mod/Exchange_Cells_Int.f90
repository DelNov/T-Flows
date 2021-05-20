!==============================================================================!
  subroutine Exchange_Cells_Int(Grid, phi)
!------------------------------------------------------------------------------!
!   Exchanges the values of a cell-based integer array between the processors. !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: phi(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: ln, c1, c2, len_s, len_r, sub
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    len_s = Grid % Comm % cells_send(sub) % n_items
    do ln = 1, len_s
      c1 = Grid % Comm % cells_send(sub) % map(ln)
      Grid % Comm % cells_send(sub) % i_buff(ln) = phi(c1)
    end do
  end do

  ! Exchange the values
  do sub = 1, n_proc
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
  do sub = 1, n_proc
    len_r = Grid % Comm % cells_recv(sub) % n_items
    do ln = 1, len_r
      c2 = Grid % Comm % cells_recv(sub) % map(ln)
      phi(c2) = Grid % Comm % cells_recv(sub) % i_buff(ln)
    end do
  end do

  end subroutine
