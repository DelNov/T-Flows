!==============================================================================!
  subroutine Control_Mod_Mass_Flow_Rates(m_x, m_y, m_z, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: m_x, m_y, m_z
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('MASS_FLOW_RATES', 3, def,  &
                                    val, verbose)

  m_x = val(1)
  m_y = val(2)
  m_z = val(3)

  end subroutine
