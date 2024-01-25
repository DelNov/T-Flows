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
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e, vapour, liquid
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  if(Flow % mass_transfer_model == 0) return

  ! Distinguish between liquid and vapor
  call Vof % Get_Vapour_And_Liquid_Phase(vapour, liquid)

  ! Add volume over all cells, avoiding buffer cells
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    e = Vof % Front % elem_in_cell(c)  ! front element
    if(e .ne. 0) then

      ! Unit: kg/s * m^3/kg = m^3/s
      b(c) = b(c)                       &
           + Vof % m_dot(c)             &
           * (1.0/Vof % phase_dens(vapour))

    end if
  end do

  end subroutine
