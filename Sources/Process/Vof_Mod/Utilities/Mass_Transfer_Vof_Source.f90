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
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e, g, l
!==============================================================================!

  ! Take aliases
  grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Distinguish between liquid and vapor
  call Vof % Get_Gas_And_Liquid_Phase(g, l)

  !-------------------!
  !   Volume source   !
  !-------------------!

  ! Here is the trick to get the sign correct:
  ! - if gas is 1 and liquid 2 => l-g =  1 => source > 0
  ! - if gas is 2 and liquid 1 => l-g = -1 => source < 0
  do c = 1, grid % n_cells
    e = Vof % Front % cell_at_elem(c)  ! Front element

    ! As Yohei and Lubomir ademantly told me - you divide with the density
    ! of the phase for "which you are solving".  And you are solving for
    ! the one which is defined as one (not zero) in the system
    if(e .ne. 0) then
      b(c) = b(c)                                                    &
           + Vof % m_dot(c) * (l-g) * Vof % Front % elem(e) % area   &
           * (1.0/Vof % phase_dens(l))
!          * (1.0/Vof % phase_dens(g) - 1.0/Vof % phase_dens(l))
    end if
  end do

  end subroutine
