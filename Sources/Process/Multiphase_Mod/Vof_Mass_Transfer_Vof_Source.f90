!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer_Vof_Source(mult, b)
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
  type(Grid_Type), pointer :: grid
  integer                  :: c, e, g, l
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid

  ! Array m_dot is positive when boiling
  mult % m_dot(:) = 0.0
  do c = 1, grid % n_cells
    if(mult % cell_at_elem(c) .ne. 0) then
      mult % m_dot(c) = 0.05  ! set it to be the velocity
    end if
  end do

  ! Distinguish between liquid and vapor
  call Multiphase_Mod_Vof_Get_Gas_And_Liquid_Phase(mult, g, l)

  ! Volume source
  ! Here is the trick to get the sign correct:
  ! - if gas is 1 and liquid 2 => l-g =  1 => source > 0
  ! - if gas is 2 and liquid 1 => l-g = -1 => source < 0
  do c = 1, grid % n_cells
    e = mult % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      b(c) = b(c)            &
           + mult % m_dot(c) * (l-g) * mult % front % elem(e) % area
    end if
  end do

  end subroutine
