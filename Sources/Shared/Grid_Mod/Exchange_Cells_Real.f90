!==============================================================================!
  subroutine Grid_Mod_Exchange_Cells_Real(grid, phi)
!------------------------------------------------------------------------------!
!   Exchanges the values of a cell-based real array between the processors.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: ln, c1, c2, len_s, len_r, sub
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    len_s = grid % comm % cells_send(sub) % n_items
    do ln = 1, len_s
      c1 = grid % comm % cells_send(sub) % map(ln)
      grid % comm % cells_send(sub) % r_buff(ln) = phi(c1)
    end do
  end do

  ! Exchange the values
  do sub = 1, n_proc
    len_s = grid % comm % cells_send(sub) % n_items
    len_r = grid % comm % cells_recv(sub) % n_items
    if(len_s + len_r > 0) then
      call Comm_Mod_Sendrecv_Real_Arrays(           &
        len_s,                                      &  ! sending length
        grid % comm % cells_send(sub) % r_buff(1),  &  ! array to be sent
        len_r,                                      &  ! receiving length
        grid % comm % cells_recv(sub) % r_buff(1),  &  ! array to be received
        sub)                                           ! destination processor
    end if
  end do

  ! Fill the buffers with new values
  do sub = 1, n_proc
    len_r = grid % comm % cells_recv(sub) % n_items
    do ln = 1, len_r
      c2 = grid % comm % cells_recv(sub) % map(ln)
      phi(c2) = grid % comm % cells_recv(sub) % r_buff(ln)
    end do
  end do

  end subroutine
