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
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! Integrate added volume
  added_vol = 0.0

  if(.not. Flow % mass_transfer) return

  ! Integrated added volume over all cells, avoiding buffer cells
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    e = Vof % Front % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      added_vol = added_vol  &
                + Vof % m_dot(c) * Vof % Front % elem(e) % area

    end if
  end do

  ! Take global summ
  call Comm_Mod_Global_Sum_Real(added_vol)

  end subroutine
