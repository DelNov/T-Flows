!==============================================================================!
  subroutine Grid_Mod_Exchange_Real(grid, phi)
!------------------------------------------------------------------------------!
!   Exchanges the values of a real array between the processors.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, sub
!==============================================================================!

  ! Fill the buffers with new values
  do sub = 1, n_proc
    if( grid % comm % nbb_e(sub)  >=   &
        grid % comm % nbb_s(sub) ) then
      do c2 = grid % comm % nbb_s(sub),  &
              grid % comm % nbb_e(sub)
        c1 = grid % comm % buffer_index(c2)
        phi(c2) = phi(c1)
      end do
    end if
  end do

  ! Exchange the values
  do sub = 1, n_proc
    if( grid % comm % nbb_e(sub)  >=   &
        grid % comm % nbb_s(sub) ) then

      call Comm_Mod_Exchange_Real_Array(       &
        phi(grid % comm % nbb_s(sub)),         &  ! array to be exchanged
          grid % comm % nbb_e(sub)             &
        - grid % comm % nbb_s(sub) + 1,        &  ! array's length
        sub)                                      ! destination processor

    end if
  end do

  end subroutine
