!==============================================================================!
  subroutine Mass_Transfer_Vof_Source(Vof, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
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

  ! If not a problem with mass transfer, get out of here
  if(Flow % mass_transfer_model == 0) return

  ! Distinguish between liquid and vapor
  call Vof % Get_Vapour_And_Liquid_Phase(vapour, liquid)

  !-------------------!
  !   Volume source   !
  !-------------------!

  ! Here is the trick to get the sign correct:
  ! - if gas is 1 and liquid 2 => liquid-vapour =  1 => source > 0
  ! - if gas is 2 and liquid 1 => liquid-vapour = -1 => source < 0
  do c = 1, Grid % n_cells
    e = Vof % Front % elem_in_cell(c)  ! Front element

    ! As Yohei and Lubomir ademantly told me - you divide with the density
    ! of the phase for "which you are solving".  And you are solving for
    ! the one which is defined as one (not zero) in the system
    !
    ! I have certain reservations about the sign (liquid-vapour).  Since the
    ! term is divided by liquid density, it is very small by definition.
    ! Also, when I changed definition of phases from 1:2 to 0:1, the results
    ! were not affected at all, but they should have been.
    if(e .ne. 0) then
      b(c) = b(c)                              &
           + Vof % m_dot(c) * (liquid-vapour)  &
           * (1.0/Vof % phase_dens(liquid))
!          * (1.0/Vof % phase_dens(vapour) - 1.0/Vof % phase_dens(liquid))
    end if
  end do

  end subroutine
