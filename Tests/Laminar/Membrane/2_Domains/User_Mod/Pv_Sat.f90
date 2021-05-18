!==============================================================================!
  subroutine Pv_Sat(t, p_v_h2o)
!------------------------------------------------------------------------------!
!   Saturation vapor pressure according to temperature and salt concentration  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(in)  :: t
  real, intent(out) :: p_v_h2o  ! partial vapor pressure in air domain
!-----------------------------------[Locals]-----------------------------------!
  real :: t_tmp  ! temperature in K
!==============================================================================!

  t_tmp = t + 273.15 ! convert to K
  if(t_tmp .gt. 334.0 .and. t_tmp .le. 363.15) then
    p_v_h2o = 1E+5 * 10**(5.0768  - 1659.793/(t_tmp-45.854))
  else if (t_tmp .gt. 304 .and. t_tmp .le. 334.0) then
    p_v_h2o = 1E+5 * 10**(5.20389 - 1733.926/(t_tmp-39.485))
  else if (t_tmp .gt. 273 .and. t_tmp .le. 304.0) then
    p_v_h2o = 1E+5 * 10**(5.40221 - 1838.675/(t_tmp-31.737))
  else
    error stop "temperature value out of saturation pressure range [273-363K]"
  end if

  end subroutine
