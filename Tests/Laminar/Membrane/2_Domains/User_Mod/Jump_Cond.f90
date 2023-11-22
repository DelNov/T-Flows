!#include "Pv_Sat.f90"

!==============================================================================!
  subroutine Jump_Cond(t_int, res, lhs_lin, lhs_fun, rhs)
!------------------------------------------------------------------------------!
!   Saturation vapor pressure according to temperature and salt concentration  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: t_int         ! interface temperature
  real :: res           ! residual (ideally zero)
  real :: lhs_lin, lhs_fun, rhs
!-----------------------------------[Locals]-----------------------------------!
  real :: pv_tint  ! partial vapor pressure
!==============================================================================!

  call Pv_Sat(t_int, pv_tint)

  res = lhs_lin * t_int + lhs_fun * pv_tint - rhs

  end subroutine
