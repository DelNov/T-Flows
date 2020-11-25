!==============================================================================!
  subroutine Grid_Mod_Exchange_Nodes_Log(grid, phi)
!------------------------------------------------------------------------------!
!   Exchanges the values of a node-based logical array between the processors. !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: phi(1:grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: ln, n, length, sub
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    length = grid % comm % nodes_repl(sub) % n_items
    do ln = 1, length
      n = grid % comm % nodes_repl(sub) % map(ln)
      grid % comm % nodes_repl(sub) % l_buff(ln) = phi(n)
    end do
  end do

  ! Exchange the values
  do sub = 1, n_proc
    length = grid % comm % nodes_repl(sub) % n_items
    if( length > 0 ) then
      call Comm_Mod_Exchange_Log_Array(             &
        length,                                     &  ! array's length
        grid % comm % nodes_repl(sub) % l_buff(1),  &  ! array to be exchanged
        sub)                                           ! destination processor
    end if
  end do

  end subroutine
