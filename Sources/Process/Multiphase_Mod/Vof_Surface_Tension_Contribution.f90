!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
!------------------------------------------------------------------------------!
!        Computes the Curvature based on old-fashion  Brackbill's CSF          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum_v1 => r_cell_01,  &
                      sum_v2 => r_cell_02,  &
                      sum_v3 => r_cell_03,  &
                      sum_v4 => r_cell_04,  &
                      sum_v5 => r_cell_05,  &
                      sum_v6 => r_cell_06
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),pointer :: flow
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: s, c, c1, c2, c_iter, n_conv
  real                     :: vol_face, grad_face(3)
  real                     :: dotprod, sxyz_mod, fs, epsloc, err_2
!==============================================================================!

  err_2 = 1.0e-14
  epsloc = epsilon(epsloc)

  ! First take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid
  vof  => mult % vof

  do c = -grid % n_bnd_cells, grid % n_cells
    mult % curv(c) = vof % n(c)
  end do

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
        mult % curv(c2) = mult % curv(c1)
      end if
    end do

    ! Smoothing vof:
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      sxyz_mod = sqrt(grid % sx(s) ** 2 + grid % sy(s) ** 2 + grid % sz(s) ** 2)
      if (c2 > 0) then
        vol_face = fs * mult % curv(c1) + (1.0 - fs) * mult % curv(c2)
        sum_v1(c1) = sum_v1(c1) + vol_face * sxyz_mod
        sum_v2(c1) = sum_v2(c1) + sxyz_mod 
        sum_v1(c2) = sum_v1(c2) + vol_face * sxyz_mod
        sum_v2(c2) = sum_v2(c2) + sxyz_mod 
      else
        vol_face = mult % curv(c1)
        sum_v1(c1) = sum_v1(c1) + vol_face * sxyz_mod
        sum_v2(c1) = sum_v2(c1) + sxyz_mod 
      end if
    end do

    call Comm_Mod_Exchange_Real(grid, sum_v1)
    call Comm_Mod_Exchange_Real(grid, sum_v2)

    do c = 1, grid % n_cells
      mult % curv(c) = sum_v1(c) / sum_v2(c)
    end do

    call Comm_Mod_Exchange_Real(grid, mult % curv)

!    if (c_iter == 6) then
!      vof % min = mult % curv
!    end if
  end do

  !clean noise
  do c = 1, grid % n_cells
    if (mult % curv(c) < epsloc) then
      mult % curv(c) = 0.0
    end if
    if ((1.0 - mult % curv(c)) < epsloc) then
      mult % curv(c) = 1.0
    end if
  end do 

  !-------------------!
  !   Find Gradient   !
  !-------------------!

  sum_v1 = 0.0
  sum_v2 = 0.0
  sum_v3 = 0.0
  sum_v4 = 0.0
  sum_v5 = 0.0
  sum_v6 = 0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (c2 > 0) then
      vol_face = fs * mult % curv(c1) + (1.0 - fs) * mult % curv(c2)

      sum_v1(c1) = sum_v1(c1) + vol_face * grid % sx(s)
      sum_v2(c1) = sum_v2(c1) + vol_face * grid % sy(s)
      sum_v3(c1) = sum_v3(c1) + vol_face * grid % sz(s)

      sum_v1(c2) = sum_v1(c2) - vol_face * grid % sx(s)
      sum_v2(c2) = sum_v2(c2) - vol_face * grid % sy(s)
      sum_v3(c2) = sum_v3(c2) - vol_face * grid % sz(s)

!      vol_face = fs * vof % min(c1) + (1.0 - fs) * vof % min(c2)
!      sum_v4(c1) = sum_v4(c1) + vol_face * grid % sx(s)
!      sum_v5(c1) = sum_v5(c1) + vol_face * grid % sy(s)
!      sum_v6(c1) = sum_v6(c1) + vol_face * grid % sz(s)
!
!      sum_v4(c2) = sum_v4(c2) - vol_face * grid % sx(s)
!      sum_v5(c2) = sum_v5(c2) - vol_face * grid % sy(s)
!      sum_v6(c2) = sum_v6(c2) - vol_face * grid % sz(s)
    else
      vol_face = mult % curv(c1)

      sum_v1(c1) = sum_v1(c1) + vol_face * grid % sx(s)
      sum_v2(c1) = sum_v2(c1) + vol_face * grid % sy(s)
      sum_v3(c1) = sum_v3(c1) + vol_face * grid % sz(s)

!      vol_face = vof % min(c1)
!      sum_v4(c1) = sum_v4(c1) + vol_face * grid % sx(s)
!      sum_v5(c1) = sum_v5(c1) + vol_face * grid % sy(s)
!      sum_v6(c1) = sum_v6(c1) + vol_face * grid % sz(s)
    end if
  end do

  do c = 1, grid % n_cells
    !clean noise
    if (abs(sum_v1(c)) < epsloc) then
      sum_v1(c) = 0.0
    end if
    if (abs(sum_v2(c)) < epsloc) then
      sum_v2(c) = 0.0
    end if
    if (abs(sum_v3(c)) < epsloc) then
      sum_v3(c) = 0.0
    end if
    sum_v1(c) = sum_v1(c) / grid % vol(c)
    sum_v2(c) = sum_v2(c) / grid % vol(c)
    sum_v3(c) = sum_v3(c) / grid % vol(c)

!    sum_v4(c) = sum_v4(c) / grid % vol(c)
!    sum_v5(c) = sum_v5(c) / grid % vol(c)
!    sum_v6(c) = sum_v6(c) / grid % vol(c)
  end do

  call Comm_Mod_Exchange_Real(grid, sum_v1)
  call Comm_Mod_Exchange_Real(grid, sum_v2)
  call Comm_Mod_Exchange_Real(grid, sum_v3)

!  call Comm_Mod_Exchange_Real(grid, sum_v4)
!  call Comm_Mod_Exchange_Real(grid, sum_v5)
!  call Comm_Mod_Exchange_Real(grid, sum_v6)

  !vof % x = sum_v4
  !vof % y = sum_v5
  !vof % z = sum_v6

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % curv = 0.0

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

      if (sxyz_mod > epsloc) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        mult % curv(c1) = mult % curv(c1) + dotprod
        mult % curv(c2) = mult % curv(c2) - dotprod
      end if

    else
      grad_face(1) = sum_v1(c1)
      grad_face(2) = sum_v2(c1)
      grad_face(3) = sum_v3(c1)

      sxyz_mod = sqrt(grad_face(1) ** 2  &
                    + grad_face(2) ** 2  &
                    + grad_face(3) ** 2)

      if (sxyz_mod > epsloc) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        mult % curv(c1) = mult % curv(c1) + dotprod
      end if

    end if
  end do

  call Comm_Mod_Exchange_Real(grid, mult % curv)


  do c = 1, grid % n_cells
    mult % curv(c) = - mult % curv(c) / grid % vol(c)
  end do

  !clean noise
  do c = 1, grid % n_cells
    if (abs(mult % curv(c)) < epsloc) then
      mult % curv(c) = 0.0
    end if
  end do 

  end subroutine
