!==============================================================================!
  subroutine Comm_Mod_Exchange_Int(grid, phi)
!------------------------------------------------------------------------------!
!   Exchanges the values of an integer array between the processors.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: phi(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, sub, rtag, stag, length, error
  integer :: status(MPI_STATUS_SIZE)
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    if( grid % comm % nbb_e(sub)  <=   &
        grid % comm % nbb_s(sub) ) then
      do c2 = grid % comm % nbb_s(sub),  &
              grid % comm % nbb_e(sub), -1
        c1 = grid % comm % buffer_index(c2)
        phi(c2) = phi(c1)
      end do
    end if
  end do

  ! Exchange the values
  do sub = 1, n_proc
    if( grid % comm % nbb_e(sub)  <=   &
        grid % comm % nbb_s(sub) ) then

      length = grid % comm % nbb_s(sub)  &
             - grid % comm % nbb_e(sub) + 1
      stag = (n_proc) * this_proc + sub  ! tag for sending
      rtag = (n_proc) * sub + this_proc  ! tag for receiving

      call Mpi_Sendrecv_Replace(phi(grid % comm % nbb_e(sub)),    & ! buffer
                                length,                           & ! length
                                MPI_INTEGER8,                     & ! datatype
                                (sub-1),                          & ! dest,
                                stag,                             & ! sendtag,
                                (sub-1),                          & ! source,
                                rtag,                             & ! recvtag,
                                MPI_COMM_WORLD,                   &
                                status,                           &
                                error)

    end if  !  nbb_e(sub) .ne. nbb_s(sub)
  end do

  end subroutine
