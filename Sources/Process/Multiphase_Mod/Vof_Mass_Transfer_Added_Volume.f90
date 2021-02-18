!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer_Added_Volume(mult, added_vol)
!------------------------------------------------------------------------------!
!   Calculates added volume due to phase change                                !
!                                                                              !
!   Called from Balance_Volume                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: added_vol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, e
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid

  ! Integrate added volume
  added_vol = 0.0
  if(mult % mass_transfer) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      e = mult % cell_at_elem(c)  ! front element
      if(e .ne. 0) then
        added_vol = added_vol  &
                  + mult % m_dot(c) * mult % front % elem(e) % area

      end if
    end do
  end if

  ! Take global summ
  call Comm_Mod_Global_Sum_Real(added_vol)

  end subroutine
