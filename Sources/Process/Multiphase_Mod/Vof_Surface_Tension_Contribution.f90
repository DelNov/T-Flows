!==============================================================================!
  subroutine Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
!------------------------------------------------------------------------------!
!    Computes the Curvature based on old-fashioned Brackbill's CSF approach    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: grad_kx   => r_cell_02,  & !grad on x of vof for curvature
                      grad_ky   => r_cell_03,  & !grad on y of vof for curvature
                      grad_kz   => r_cell_04,  & !grad on z of vof for curvature
                      grad_nx   => r_cell_05,  & !grad on x of vof for normal
                      grad_ny   => r_cell_06,  & !grad on y of vof for normal
                      grad_nz   => r_cell_07,  & !grad on z of vof for normal
                      vof_norm  => r_cell_08   ! normal
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),pointer :: flow
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: vof
  integer                  :: s, c, c1, c2, c_iter
  real                     :: vol_face, grad_face(3), grad_control(3)
  real                     :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                     :: d_n(3) !normal pointing to the wall
  real                     :: norm_grad !normal of a gradient
!==============================================================================!

  epsloc = epsilon(epsloc)

  ! First take aliases
  flow => mult % pnt_flow
  grid => mult % pnt_grid
  vof  => mult % vof

  if(mult % d_func) then  ! using distance function

    grad_kx = mult % dist_func % x
    grad_ky = mult % dist_func % y
    grad_kz = mult % dist_func % z

  else ! using VOF

    do c = -grid % n_bnd_cells, grid % n_cells
      mult % curv(c) = vof % n(c)
    end do

    do c_iter = 1, mult % n_conv_curv

      grad_kx = 0.0
      grad_ky = 0.0

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

      ! Smoothing VOF
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        fs = grid % f(s)
        sxyz_mod = sqrt( grid % sx(s) ** 2   &
                       + grid % sy(s) ** 2   &
                       + grid % sz(s) ** 2)
        if (c2 > 0) then
          vol_face = fs * mult % curv(c1) + (1.0 - fs) * mult % curv(c2)
          grad_kx(c1) = grad_kx(c1) + vol_face * sxyz_mod
          grad_ky(c1) = grad_ky(c1) + sxyz_mod
          grad_kx(c2) = grad_kx(c2) + vol_face * sxyz_mod
          grad_ky(c2) = grad_ky(c2) + sxyz_mod
        else
          vol_face = mult % curv(c1)
          grad_kx(c1) = grad_kx(c1) + vol_face * sxyz_mod
          grad_ky(c1) = grad_ky(c1) + sxyz_mod
        end if
      end do

      call Grid_Mod_Exchange_Real(grid, grad_kx)
      call Grid_Mod_Exchange_Real(grid, grad_ky)

      do c = 1, grid % n_cells
        mult % curv(c) = grad_kx(c) / grad_ky(c)
      end do

      call Grid_Mod_Exchange_Real(grid, mult % curv)

      if (c_iter == mult % n_conv_norm) then  ! smooth for normal
        vof_norm = mult % curv
      end if
    end do

    ! Clean noise
    do c = 1, grid % n_cells
      if (mult % curv(c) < epsloc) then
        mult % curv(c) = 0.0
      end if
      if ((1.0 - mult % curv(c)) < epsloc) then
        mult % curv(c) = 1.0
      end if
    end do 

    ! At boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if (c2 < 0) then
        mult % curv(c2) = mult % curv(c1)
        vof_norm(c2) = vof_norm(c1)
      end if
    end do

    !-------------------!
    !   Find Gradient   !
    !-------------------!

    grad_kx = 0.0
    grad_ky = 0.0
    grad_kz = 0.0
    grad_nx = 0.0
    grad_ny = 0.0
    grad_nz = 0.0

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      if (c2 > 0) then
        vol_face = fs * mult % curv(c1) + (1.0 - fs) * mult % curv(c2)
        grad_kx(c1) = grad_kx(c1) + vol_face * grid % sx(s)
        grad_ky(c1) = grad_ky(c1) + vol_face * grid % sy(s)
        grad_kz(c1) = grad_kz(c1) + vol_face * grid % sz(s)

        grad_kx(c2) = grad_kx(c2) - vol_face * grid % sx(s)
        grad_ky(c2) = grad_ky(c2) - vol_face * grid % sy(s)
        grad_kz(c2) = grad_kz(c2) - vol_face * grid % sz(s)

        vol_face = fs * vof_norm(c1) + (1.0 - fs) * vof_norm(c2)
        grad_nx(c1) = grad_nx(c1) + vol_face * grid % sx(s)
        grad_ny(c1) = grad_ny(c1) + vol_face * grid % sy(s)
        grad_nz(c1) = grad_nz(c1) + vol_face * grid % sz(s)

        grad_nx(c2) = grad_nx(c2) - vol_face * grid % sx(s)
        grad_ny(c2) = grad_ny(c2) - vol_face * grid % sy(s)
        grad_nz(c2) = grad_nz(c2) - vol_face * grid % sz(s)
      else
        vol_face = mult % curv(c1)
        grad_kx(c1) = grad_kx(c1) + vol_face * grid % sx(s)
        grad_ky(c1) = grad_ky(c1) + vol_face * grid % sy(s)
        grad_kz(c1) = grad_kz(c1) + vol_face * grid % sz(s)

        vol_face = vof_norm(c1)
        grad_nx(c1) = grad_nx(c1) + vol_face * grid % sx(s)
        grad_ny(c1) = grad_ny(c1) + vol_face * grid % sy(s)
        grad_nz(c1) = grad_nz(c1) + vol_face * grid % sz(s)
      end if
    end do

    do c = 1, grid % n_cells
      grad_kx(c) = grad_kx(c) / grid % vol(c)
      grad_ky(c) = grad_ky(c) / grid % vol(c)
      grad_kz(c) = grad_kz(c) / grid % vol(c)

      grad_nx(c) = grad_nx(c) / grid % vol(c)
      grad_ny(c) = grad_ny(c) / grid % vol(c)
      grad_nz(c) = grad_nz(c) / grid % vol(c)
    end do

  end if

  ! Imposing contact angle at walls (before finding curvature):

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    if (c2 < 0) then
      if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then !contact angle
        norm_grad = norm2((/ grad_kx(c1), grad_ky(c1), grad_kz(c1)/))
        if (norm_grad > epsloc) then
          d_n = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/)    &
                           ,(/grid % sx(s), grid % sy(s), grid % sz(s)/))   &
                           / grid % s(s) ** 2                               &
                           * (/grid % sx(s), grid % sy(s), grid % sz(s)/)

          grad_kx(c1) = grid % dx(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_kx(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad
          grad_ky(c1) = grid % dy(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_ky(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad
          grad_kz(c1) = grid % dz(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_kz(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad

        end if

        norm_grad = norm2((/ grad_nx(c1), grad_ny(c1), grad_nz(c1)/))

        if (norm_grad > epsloc) then
          d_n = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/)    &
                           ,(/grid % sx(s), grid % sy(s), grid % sz(s)/))   &
                           / grid % s(s) ** 2                               &
                           * (/grid % sx(s), grid % sy(s), grid % sz(s)/)

          grad_nx(c1) = grid % dx(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_nx(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad
          grad_ny(c1) = grid % dy(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_ny(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad
          grad_nz(c1) = grid % dz(s) / norm2(d_n)                      &
                                     * cos(vof % q(c2) * PI / 180.0)   &
                       + grad_nz(c1) * sin(vof % q(c2) * PI / 180.0) / norm_grad

        end if

      end if

    end if
  end do

  call Grid_Mod_Exchange_Real(grid, grad_kx)
  call Grid_Mod_Exchange_Real(grid, grad_ky)
  call Grid_Mod_Exchange_Real(grid, grad_kz)

  call Grid_Mod_Exchange_Real(grid, grad_nx)
  call Grid_Mod_Exchange_Real(grid, grad_ny)
  call Grid_Mod_Exchange_Real(grid, grad_nz)

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % curv = 0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (c2 > 0) then
      grad_face(1) = fs * grad_kx(c1) + (1.0 - fs) * grad_kx(c2)
      grad_face(2) = fs * grad_ky(c1) + (1.0 - fs) * grad_ky(c2)
      grad_face(3) = fs * grad_kz(c1) + (1.0 - fs) * grad_kz(c2)

      sxyz_mod = sqrt(grad_face(1) ** 2 + grad_face(2) ** 2 + grad_face(3) ** 2)

      if (sxyz_mod > epsloc) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        sxyz_control = norm2((/grad_kx(c1), grad_ky(c1), grad_kz(c1)/))

        if (sxyz_control > epsloc) then
          mult % curv(c1) = mult % curv(c1) + dotprod
        end if

        sxyz_control = norm2((/grad_kx(c2), grad_ky(c2), grad_kz(c2)/))

        if (sxyz_control > epsloc) then
          mult % curv(c2) = mult % curv(c2) - dotprod
        end if
      end if

    else
      grad_face(1) = grad_kx(c1)
      grad_face(2) = grad_ky(c1)
      grad_face(3) = grad_kz(c1)

      sxyz_mod = sqrt(grad_face(1) ** 2 + grad_face(2) ** 2 + grad_face(3) ** 2)

      if (sxyz_mod > epsloc) then
        dotprod = (grad_face(1) * grid % sx(s)  &
                 + grad_face(2) * grid % sy(s)  &
                 + grad_face(3) * grid % sz(s)) / sxyz_mod

        sxyz_control = norm2((/grad_kx(c1), grad_ky(c1), grad_kz(c1)/))

        if (sxyz_control > epsloc) then
          mult % curv(c1) = mult % curv(c1) + dotprod
        end if
      end if
    end if
  end do

  call Grid_Mod_Exchange_Real(grid, mult % curv)

  do c = 1, grid % n_cells
    mult % curv(c) = - mult % curv(c) / grid % vol(c)
  end do

  ! Smoothed normal
  if(.not.(mult % d_func) .and. mult % n_conv_norm > 0   &
     .and. mult % n_conv_norm <= mult % n_conv_curv) then
    vof % x = grad_nx
    vof % y = grad_ny
    vof % z = grad_nz
  end if

  end subroutine
