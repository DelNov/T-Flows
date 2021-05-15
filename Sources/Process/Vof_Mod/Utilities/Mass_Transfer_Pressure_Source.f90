!==============================================================================!
  subroutine Mass_Transfer_Pressure_Source(Vof, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: b(Vof % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: c, e
!==============================================================================!

  ! Take aliases
  grid => Vof % pnt_grid
  flow => Vof % pnt_flow

  RETURN

  if(.not. flow % mass_transfer) return

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    e = Vof % Front % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      b(c) = b(c)            &
           + Vof % m_dot(c) * Vof % Front % elem(e) % area
    end if
  end do

  end subroutine
