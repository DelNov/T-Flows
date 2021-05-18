!==============================================================================!
  subroutine T_Sat(t, p_v)
!------------------------------------------------------------------------------!
!   Saturation temperature from vapor pressure !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: t
  real :: p_v  ! partial vapor pressure in air domain
!-----------------------------------[Locals]-----------------------------------!
  real :: t_tmp  ! temperature in K
!==============================================================================!

  if(p_v .gt. 2.0725e+04 .and. p_v .le. 7.010428E+04) then
    t_tmp = 45.854 + 1659.793 /(5.0768  - log10(p_v * 1E-5))
  else if (p_v.gt. 4.4556e+03 .and. p_v .le. 2.0725e+04) then
    t_tmp = 39.485 + 1733.926 /(5.20389 - log10(p_v * 1E-5))
  else if (p_v .gt. 604.1854 .and. p_v .le. 4.4556e+03) then
    t_tmp = 31.737 + 1838.675 /(5.40221 - log10(p_v * 1E-5))
  else
    error stop 'saturation pressure value out of ' //  &
               'temperature calculation range [605-70''000Pa]'
  end if

  t = t_tmp - 273.15 ! convert to Celsius

  end subroutine
