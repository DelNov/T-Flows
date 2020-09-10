!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Csf(mult,  &
                                              grad_kx, grad_ky, grad_kz)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Gauss theorem        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: div_x => r_cell_10,  &
                      div_y => r_cell_11,  &
                      div_z => r_cell_12
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: grad_kx(-mult % pnt_grid % n_bnd_cells    &
                                          : mult % pnt_grid % n_cells),      &
                                   grad_ky(-mult % pnt_grid % n_bnd_cells    &
                                          : mult % pnt_grid % n_cells),      &
                                   grad_kz(-mult % pnt_grid % n_bnd_cells    &
                                          : mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: smooth
  integer                   :: s, c, c1, c2, n, i_fac,i_nod, tot_cells,sub
  integer                   :: c_inte, fu, nb, nc
  real, contiguous, pointer :: fs_x(:), fs_y(:), fs_z(:)
  real                      :: vol_face, grad_face(3), d_n(3)
  real                      :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                      :: dotprod2, stabilize
  real                      :: n_0(3), n_f(3), n_w(3), reflex(3)
  real                      :: theta, theta_0, a, b, s_vector(3)
  real                      :: vof_fx, vof_fy, vof_fz, vof_c1, vof_c2, voff
  real                      :: res1, res2, resul, term_c, sumtot
  real                      :: sumx, sumy, sumz, norm_grad, coeff
  real                      :: v1(3), v2(3), v3(3), v4(3)
  real                      :: c_c
!==============================================================================!

  grid   => mult % pnt_grid
  flow   => mult % pnt_flow
  vof    => mult % vof
  smooth => mult % smooth

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  epsloc = epsilon(epsloc)

  mult % curv = 0.0

  ! Normalize vector at cells
  do c = 1, grid % n_cells
    norm_grad = sqrt(grad_kx(c) ** 2 + grad_ky(c) ** 2 + grad_kz(c) ** 2)
    if (norm_grad >= epsloc) then
      grad_kx(c) = grad_kx(c) / norm_grad
      grad_ky(c) = grad_ky(c) / norm_grad
      grad_kz(c) = grad_kz(c) / norm_grad
    else
      grad_kx(c) = 0.0
      grad_ky(c) = 0.0
      grad_kz(c) = 0.0
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz(-nb:nc))

  ! Tangent vector to walls/symmetries

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.   &
        Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then

      norm_grad = norm2((/grad_kx(c1),grad_ky(c1),grad_kz(c1)/))

      if (norm_grad > epsloc) then
        v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
        v2 = (/grad_kx(c1), grad_ky(c1), grad_kz(c1)/)
        v3 = Math_Mod_Cross_Product(v1, v2)
        v4 = Math_Mod_Cross_Product(v3, v1)
        ! projection on v4
        norm_grad = norm2(v4)
        if (norm_grad > epsloc) then
          grad_kx(c2) = v4(1) / norm_grad
          grad_ky(c2) = v4(2) / norm_grad
          grad_kz(c2) = v4(3) / norm_grad
        end if
      end if
    else
      grad_kx(c2) = grad_kx(c1)
      grad_ky(c2) = grad_ky(c1)
      grad_kz(c2) = grad_kz(c1)
    end if
  end do

  ! Correct for contact angle at walls

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then

      ! Accumulate values of faces
      norm_grad = norm2((/grad_kx(c1),grad_ky(c1),grad_kz(c1)/))
      if (norm_grad > epsloc) then
        dotprod = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/),   &
                              (/grid % sx(s), grid % sy(s), grid % sz(s)/))

        grad_kx(c1) = grid % dx(s) / dotprod * grid % s(s)                    &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + grad_kx(c2) * sin(vof % q(c2) * PI /180.0)

        grad_ky(c1) = grid % dy(s) / dotprod * grid % s(s)                    &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + grad_ky(c2) * sin(vof % q(c2) * PI /180.0)

        grad_kz(c1) = grid % dz(s) / dotprod * grid % s(s)                    &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + grad_kz(c2) * sin(vof % q(c2) * PI /180.0)

        grad_kx(c2) = grad_kx(c1)
        grad_ky(c2) = grad_ky(c1)
        grad_kz(c2) = grad_kz(c1)
      end if

    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz(-nb:nc))

  mult % fc_x(-nb:nc) = grad_kx(-nb:nc)
  mult % fc_y(-nb:nc) = grad_ky(-nb:nc)
  mult % fc_z(-nb:nc) = grad_kz(-nb:nc)

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  ! Find divergence of normals
  call Field_Mod_Grad_Component(flow, grad_kx(-nb:nc),  &
                                1,    div_x  (-nb:nc),  &
                                impose_symmetry = .false.)
  call Field_Mod_Grad_Component(flow, grad_ky(-nb:nc),  &
                                2,    div_y  (-nb:nc),  &
                                impose_symmetry = .false.)
  call Field_Mod_Grad_Component(flow, grad_kz(-nb:nc),  &
                                3,    div_z  (-nb:nc),  &
                                impose_symmetry = .false.)

  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_x(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_y(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_z(-nb:nc)

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  call Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,                  &
                          grad_kx(-nb:nc), grad_ky(-nb:nc), grad_kz(-nb:nc))

  end subroutine
