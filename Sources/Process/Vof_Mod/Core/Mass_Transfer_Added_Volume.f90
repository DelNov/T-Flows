!==============================================================================!
  subroutine Mass_Transfer_Added_Volume(Vof, added_vol)
!------------------------------------------------------------------------------!
!   Calculates added volume due to phase change                                !
!                                                                              !
!   Called from Balance_Volume                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: added_vol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: c, e
!==============================================================================!

  ! Take aliases
  grid => Vof % pnt_grid
  flow => Vof % pnt_flow

  ! Integrate added volume
  added_vol = 0.0

  RETURN

  if(.not. flow % mass_transfer) return

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    e = Vof % Front % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      added_vol = added_vol  &
                + Vof % m_dot(c) * Vof % Front % elem(e) % area

    end if
  end do

  ! Take global summ
  call Comm_Mod_Global_Sum_Real(added_vol)

  end subroutine
