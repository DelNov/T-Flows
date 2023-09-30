# ifdef __INTEL_COMPILER
#   include "User_Mod/Pv_Sat.f90"
# else
#   include "Pv_Sat.f90"
# endif

!==============================================================================!
  subroutine Pv_Sat_Salt(t, m_h2o, m_salt, sc1, p_v_h2o)
!------------------------------------------------------------------------------!
!   Saturation vapor pressure according to temperature and salt concentration  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(in)  :: t, m_h2o, m_salt, sc1
  real, intent(out) :: p_v_h2o       ! partial vapor pressure in air domain
!-----------------------------------[Locals]-----------------------------------!
  real :: pv
!==============================================================================!

  call Pv_Sat(t, pv)
  p_v_h2o = pv * (1 - m_h2o / m_salt * sc1)

  end subroutine
