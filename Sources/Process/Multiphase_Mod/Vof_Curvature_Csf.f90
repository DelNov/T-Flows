!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Csf(grid, mult,                  &
                                              grad_kx, grad_ky, grad_kz,   &
                                              curr_colour)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Gauss theorem        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: div_x        => r_cell_10,  &
                      div_y        => r_cell_11,  &
                      div_z        => r_cell_12,  &
                      k_limit      => r_cell_13,  &
                      kappa_f      => r_face_01,  &
                      counter_nod  => i_node_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: grad_kx    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_ky    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_kz    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   curr_colour(-grid % n_bnd_cells    &
                                              : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  integer                       :: s, c, c1, c2, n, i_fac,i_nod, tot_cells,sub
  integer                       :: c_inte, fu
  real, contiguous,     pointer :: fs_x(:), fs_y(:), fs_z(:)
  real                          :: vol_face, grad_face(3), d_n(3)
  real                          :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                          :: dotprod2, stabilize
  real                          :: n_0(3), n_f(3), n_w(3), reflex(3)
  real                          :: theta, theta_0, a, b, s_vector(3)
  real                          :: vof_fx, vof_fy, vof_fz, vof_c1, vof_c2, voff
  real                          :: res1, res2, resul, term_c, sumtot
  real                          :: sumx, sumy, sumz, norm_grad, coeff
  real                          :: v1(3), v2(3), v3(3), v4(3)
  real                          :: c_c
!==============================================================================!

  vof  => mult % vof
  flow => mult % pnt_flow

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

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz)

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

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz)

  mult % fc_x = grad_kx
  mult % fc_y = grad_ky
  mult % fc_z = grad_kz

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  ! find divergence of normals

  call Multiphase_Mod_Vof_Grad_Component(flow, grad_kx, 1, div_x)
  call Multiphase_Mod_Vof_Grad_Component(flow, grad_ky, 2, div_y)
  call Multiphase_Mod_Vof_Grad_Component(flow, grad_kz, 3, div_z)

  mult % curv = mult % curv - div_x
  mult % curv = mult % curv - div_y
  mult % curv = mult % curv - div_z

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  call Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,                  &
                                           grad_kx, grad_ky, grad_kz)

  end subroutine
