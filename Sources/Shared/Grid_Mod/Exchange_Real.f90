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
    if( grid % comm % buff_e_cell(sub)  >=   &
        grid % comm % buff_s_cell(sub) ) then
      do c2 = grid % comm % buff_s_cell(sub),  &
              grid % comm % buff_e_cell(sub)
        c1 = grid % comm % buff_c1_from_c2(c2)
        phi(c2) = phi(c1)
      end do
    end if
  end do

  ! Exchange the values
  do sub = 1, n_proc
    if( grid % comm % buff_e_cell(sub)  >=   &
        grid % comm % buff_s_cell(sub) ) then

      call Comm_Mod_Exchange_Real_Array(       &
        phi(grid % comm % buff_s_cell(sub)),   &  ! array to be exchanged
          grid % comm % buff_e_cell(sub)       &  ! end minus start plus 1 ...
        - grid % comm % buff_s_cell(sub) + 1,  &  ! ... makes array's length
        sub)                                      ! destination processor

    end if
  end do

  end subroutine
