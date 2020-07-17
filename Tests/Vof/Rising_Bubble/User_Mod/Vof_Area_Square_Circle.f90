  subroutine Vof_Area_Square_Circle(xint0, xint1, cx, cy, r, h, area)
!------------------------------------------------------------------------------!
!   Computes intersection area between square and circl                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real                          :: xint0, xint1, cx, cy, r, h, area
!-----------------------------------[Locals]-----------------------------------!
  real                          :: x1, x2, xmin, xmax, tmin, tmax, xaux
  real                          :: integ1, integ2, epsloc
!==============================================================================!

  epsloc = epsilon(epsloc)

  if (xint0 > xint1) then
    xaux = xint1
    xint1 = xint0
    xint0 = xaux
  end if

  if (h > r) then
    area = 0.0
  else
    x1 = cx - sqrt(r ** 2 - h ** 2)
    x2 = cx + sqrt(r ** 2 - h ** 2)

    xmin = max(xint0, x1)
    xmax = min(xint1, x2)

    if (xmin > xmax) then
      area = 0.0
    else
      tmin = r ** 2 - (cx - xmin) ** 2 + epsloc
      tmax = r ** 2 - (cx - xmax) ** 2 + epsloc

      integ1 = - 0.5 * (cx - xmax) * sqrt(tmax)    &
               + 0.5 * r ** 2 * atan( (xmax - cx) * sqrt(tmax) / tmax) - h * xmax

      integ2 = - 0.5 * (cx - xmin) * sqrt(tmin)    &
               + 0.5 * r ** 2 * atan( (xmin - cx) * sqrt(tmin) / tmin) - h * xmin

      area = integ1 - integ2
    end if
  end if
  end subroutine
