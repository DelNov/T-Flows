!==============================================================================!
  subroutine Exchange_Cells_Real(Grid, phi, caller)
!------------------------------------------------------------------------------!
!   Exchanges the values of a cell-based real array between the processors.    !
!                                                                              !
!   In 2021, a neat extension was introduced to this procedure.  One can       !
!   send an optional parameter, a string, with name and line of the caller     !
!   procedure.  (For example: 'Compute_Momentum_347'.)  If that parameter      !
!   is sent, this procedure will check if the call was needed or reduntant.    !
!                                                                              !
!   The info is stored in file 'exchange_cells_real.log'.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)           :: Grid
  real                       :: phi(-Grid % n_bnd_cells:Grid % n_cells)
  character(len=*), optional :: caller
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

    if(this_proc < 2) then
      call File_Mod_Append_File_For_Writing(  &
                'exchange_cells_real.log', fu, this_proc)
    end if

    do sub = 1, n_proc
      len_r = Grid % comm % cells_recv(sub) % n_items
      do ln = 1, len_r
        c2 = Grid % comm % cells_recv(sub) % map(ln)
        Grid % comm % cells_recv(sub) % o_buff(ln) = phi(c2)
      end do
    end do
  end if

  !----------------------------------!
  !   Usual course of this routine   !
  !----------------------------------!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    len_s = Grid % comm % cells_send(sub) % n_items
    do ln = 1, len_s
      c1 = Grid % comm % cells_send(sub) % map(ln)
      Grid % comm % cells_send(sub) % r_buff(ln) = phi(c1)
    end do
  end do

  ! Exchange the values
  do sub = 1, n_proc
    len_s = Grid % comm % cells_send(sub) % n_items
    len_r = Grid % comm % cells_recv(sub) % n_items
    if(len_s + len_r > 0) then
      call Comm_Mod_Sendrecv_Real_Arrays(           &
        len_s,                                      &  ! sending length
        Grid % comm % cells_send(sub) % r_buff(1),  &  ! array to be sent
        len_r,                                      &  ! receiving length
        Grid % comm % cells_recv(sub) % r_buff(1),  &  ! array to be received
        sub)                                           ! destination processor
    end if
  end do

  ! Fill the buffers with new values
  do sub = 1, n_proc
    len_r = Grid % comm % cells_recv(sub) % n_items
    do ln = 1, len_r
      c2 = Grid % comm % cells_recv(sub) % map(ln)
      phi(c2) = Grid % comm % cells_recv(sub) % r_buff(ln)
    end do
  end do

  !-------------------------------------------------!
  !   Optional parameter was sent; check new with   !
  !    old values to see if this call was needed    !
  !-------------------------------------------------!
  if(present(caller)) then

    needed = NO
    do sub = 1, n_proc
      len_r = Grid % comm % cells_recv(sub) % n_items
      do ln = 1, len_r
        if( Grid % comm % cells_recv(sub) % r_buff(ln) .ne. &
            Grid % comm % cells_recv(sub) % o_buff(ln) ) then
          needed = YES
        end if
      end do
    end do

    ! What if it is needed on some other processor
    call Comm_Mod_Global_Max_Int(needed)

    if(this_proc < 2) then
      if(needed .eq. YES) then
        write(fu,'(a)') '# Call from: ' // caller // ' was needed!'
      else
        write(fu,'(a)') '# Call from: ' // caller // ' was redundant!'
      end if
      close(fu)
    end if

  end if

  end subroutine
