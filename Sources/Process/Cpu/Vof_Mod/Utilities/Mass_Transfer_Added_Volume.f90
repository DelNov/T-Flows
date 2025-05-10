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
  integer                   :: c, e, vapour, liquid
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! Integrate added volume
  added_vol = 0.0

  if(Flow % mass_transfer_model .eq. NO_MASS_TRANSFER) return

  ! Distinguish between vapour and liquid
  call Vof % Get_Vapour_And_Liquid_Phase(vapour, liquid)

  ! Integrated added volume over all cells, avoiding buffer cells
  do c = Cells_In_Domain()
    e = Vof % Front % elem_in_cell(c)  ! front element
    if(e .ne. 0) then

      ! Unit: kg/s * m^3/kg = m^3/s
      added_vol = added_vol                  &
                + Vof % m_dot(c)             &
                * (1.0/Vof % phase_dens(vapour))

    end if
  end do

  ! Take global summ
  call Global % Sum_Real(added_vol)

  end subroutine
