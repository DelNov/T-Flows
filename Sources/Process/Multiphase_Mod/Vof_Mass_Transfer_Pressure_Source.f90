!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer_Pressure_Source(mult, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: b(mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: c, e
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid
  flow => mult % pnt_flow

  RETURN

  if(.not. flow % mass_transfer) return

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    e = mult % Front % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      b(c) = b(c)            &
           + mult % m_dot(c) * mult % Front % elem(e) % area
    end if
  end do

  end subroutine
