!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
!------------------------------------------------------------------------------!
!   Computes the value at the cell face using different convective  schemes.   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_v1 => r_cell_01,  &
                      sum_v2 => r_cell_02,  &
                      sum_v3 => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, c_iter,n_conv
  real                     :: vol_face, grad_face(3)
  real                     :: dotprod, sxyz_mod, fs
!==============================================================================!

  ! First take aliases
  grid => mult % pnt_grid

  mult % vof % oo = mult % vof % n

  ! This parameter maybe also implemented in "control" file
  n_conv = 2

  do c_iter = 1, n_conv

    sum_v1 = 0.0
    sum_v2 = 0.0

    !-------------------------------!
    !   Extrapolate to boundaries   !
    !-------------------------------!

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if (c2 < 0) then
        mult % vof % oo(c2) = mult % vof % oo(c1) 
      end if
    end do

    ! Smoothing vof:
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      sxyz_mod = sqrt(grid % sx(s) ** 2 + grid % sy(s) ** 2 + grid % sz(s) ** 2)
      if (c2 > 0) then
        vol_face = fs * mult % vof % oo(c1) + (1.0 - fs) * mult % vof % oo(c2)
        sum_v1(c1) = sum_v1(c1) + vol_face * sxyz_mod
        sum_v2(c1) = sum_v2(c1) + sxyz_mod 
        sum_v1(c2) = sum_v1(c2) + vol_face * sxyz_mod
        sum_v2(c2) = sum_v2(c2) + sxyz_mod 
      else
        vol_face = mult % vof % oo(c1)
        sum_v1(c1) = sum_v1(c1) + vol_face * sxyz_mod
        sum_v2(c1) = sum_v2(c1) + sxyz_mod 
      end if
    end do

    do c = 1, grid % n_cells
      mult % vof % oo(c) = sum_v1(c) / sum_v2(c)
    end do

  end do

  !-------------------!
  !   Find Gradient   !
  !-------------------!

  sum_v1 = 0.0
  sum_v2 = 0.0
  sum_v3 = 0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (c2 > 0) then
      vol_face = fs * mult % vof % oo(c1) + (1.0 - fs) * mult % vof % oo(c2)
      sum_v1(c1) = sum_v1(c1) + vol_face * grid % sx(s)
      sum_v2(c1) = sum_v2(c1) + vol_face * grid % sy(s)
      sum_v3(c1) = sum_v3(c1) + vol_face * grid % sz(s)

      sum_v1(c2) = sum_v1(c2) - vol_face * grid % sx(s)
      sum_v2(c2) = sum_v2(c2) - vol_face * grid % sy(s)
      sum_v3(c2) = sum_v3(c2) - vol_face * grid % sz(s)
    else
      vol_face = mult % vof % oo(c1)
      sum_v1(c1) = sum_v1(c1) + vol_face * grid % sx(s)
      sum_v2(c1) = sum_v2(c1) + vol_face * grid % sy(s)
      sum_v3(c1) = sum_v3(c1) + vol_face * grid % sz(s)
    end if
  end do

  do c = 1, grid % n_cells
    sum_v1(c) = sum_v1(c) / grid % vol(c)
    sum_v2(c) = sum_v2(c) / grid % vol(c)
    sum_v3(c) = sum_v3(c) / grid % vol(c)
  end do

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % vof % oo =0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (c2 > 0) then
      grad_face(1) = fs * sum_v1(c1) + (1.0 - fs) * sum_v1(c2)
      grad_face(2) = fs * sum_v2(c1) + (1.0 - fs) * sum_v2(c2)
      grad_face(3) = fs * sum_v3(c1) + (1.0 - fs) * sum_v3(c2)

      sxyz_mod = sqrt(grad_face(1) ** 2  &
                    + grad_face(2) ** 2  &
                    + grad_face(3) ** 2)

      if (sxyz_mod > TINY) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        mult % vof % oo(c1) = mult % vof % oo(c1) + dotprod
        mult % vof % oo(c2) = mult % vof % oo(c2) - dotprod
      else
        mult % vof % oo(c1) = 0.0
        mult % vof % oo(c2) = 0.0
      end if

    else
      grad_face(1) = sum_v1(c1)
      grad_face(2) = sum_v2(c1)
      grad_face(3) = sum_v3(c1)

      sxyz_mod = sqrt(grad_face(1) ** 2  &
                    + grad_face(2) ** 2                   &
                    + grad_face(3) ** 2)

      if (sxyz_mod > TINY) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        mult % vof % oo(c1) = mult % vof % oo(c1) + dotprod
      else
        mult % vof % oo(c1) = 0.0
      end if

    end if
  end do

  do c = 1, grid % n_cells
    mult % vof % oo(c) = - mult % vof % oo(c) / grid % vol(c)
  end do

  end subroutine
