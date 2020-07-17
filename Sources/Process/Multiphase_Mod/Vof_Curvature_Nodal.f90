!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Nodal(grid, mult, smooth_k,   &
                                                var_node_k)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using nodal gradients      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: div_x  => r_cell_11,  &
                      div_y  => r_cell_12,  &
                      div_z  => r_cell_13,  &
                      div_xx => r_cell_14,  &
                      div_yy => r_cell_15,  &
                      div_zz => r_cell_16,  &
                      grad_x => r_node_05,  &
                      grad_y => r_node_06,  &
                      grad_z => r_node_07,  &
                      mark   => i_node_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: smooth_k(-grid % n_bnd_cells :   &
                                             grid % n_cells),       &
                                   var_node_k(1:grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  integer                       :: s, c, c1, c2, n, i_fac, i_nod, tot_cells,sub
  integer                       :: c_inte, fu, nb, nc, nn
  real, contiguous,     pointer :: fs_x(:), fs_y(:), fs_z(:)
  real                          :: vol_face, grad_face(3), d_n(3)
  real                          :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                          :: dotprod2, stabilize
  real                          :: n_0(3), n_f(3), n_w(3), reflex(3)
  real                          :: theta, theta0, a, b, s_vector(3)
  real                          :: vof_fx, vof_fy, vof_fz, vof_c1, vof_c2, voff
  real                          :: res1, res2, resul, term_c, sumtot
  real                          :: sumx, sumy, sumz, norm_grad, coeff
  real                          :: v1(3), v2(3), v3(3), v4(3)
  real                          :: c_c
!==============================================================================!

  vof  => mult % vof
  flow => mult % pnt_flow

  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes

  epsloc = epsilon(epsloc)

  ! Calculate gradients at nodes
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            1, grad_x(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            2, grad_y(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            3, grad_z(1:nn))

  ! Tangent vector to symmetries
  mark(1:nn) = 0

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(     (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY)    &
       .or. (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) ) then

      do i_nod = 1, grid % faces_n_nodes(s)
        n = grid % faces_n(i_nod,s)
        if(mark(n) == 0) then
          v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
          v2 = (/grad_x(n), grad_y(n), grad_z(n)/)
          v3 = Math_Mod_Cross_Product(v1, v2)
          v4 = Math_Mod_Cross_Product(v3, v1)
          ! projection on v4
          norm_grad = norm2(v4)
          if (norm_grad > epsloc) then
            grad_x(n) = v4(1) / norm_grad
            grad_y(n) = v4(2) / norm_grad
            grad_z(n) = v4(3) / norm_grad
            mark(n) = 1
          end if
        end if
      end do

    end if
  end do

  ! Normalize vectors at nodes
  do n = 1, grid % n_nodes
    norm_grad = sqrt(grad_x(n) ** 2 + grad_y(n) ** 2 + grad_z(n) ** 2)
    if (norm_grad >= epsloc) then
      grad_x(n) = grad_x(n) / norm_grad
      grad_y(n) = grad_y(n) / norm_grad
      grad_z(n) = grad_z(n) / norm_grad
    else
      grad_x(n) = 0.0
      grad_y(n) = 0.0
      grad_z(n) = 0.0
    end if
  end do

  ! Contact angle
  mark(1:nn) = 0

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then

      do i_nod = 1, grid % faces_n_nodes(s)
        n = grid % faces_n(i_nod,s)
        if(mark(n) == 0) then
          v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/) / grid % s(s)
          v2 = (/grad_x(n), grad_y(n), grad_z(n)/)
          theta = vof % q(c2) * PI /180.0
          theta0 = acos(dot_product(v1,v2))

          a = (cos(theta) - cos(theta0) * cos(theta0 - theta))    &
            / (1.0 - cos(theta0) ** 2)
          b = (cos(theta0 - theta) - cos(theta0) * cos(theta))    &
            / (1.0 - cos(theta0) ** 2)

          norm_grad = norm2(v2)
          if (norm_grad > epsloc) then
            grad_x(n) = v2(1) * a + v1(1) * b
            grad_y(n) = v2(2) * a + v1(2) * b
            grad_z(n) = v2(3) * a + v1(3) * b
            mark(n) = 1
          end if
        end if
      end do

    end if
  end do

  call Grid_Mod_Exchange_Nodes_Real(grid, grad_x(1:nn))
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_y(1:nn))
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_z(1:nn))

  ! Interpolate node values to cells
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_x(1:nn), div_x(-nb:nc))
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_y(1:nn), div_y(-nb:nc))
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_z(1:nn), div_z(-nb:nc))

  ! Correct for contact angle at walls
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then

      ! Accumulate values of faces
      norm_grad = norm2((/div_x(c1),div_y(c1),div_z(c1)/))

      if (norm_grad > epsloc) then
        dotprod = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/),   &
                              (/grid % sx(s), grid % sy(s), grid % sz(s)/))

        div_x(c1) = grid % dx(s) / dotprod * grid % s(s)         &
                    * cos(vof % q(c2) * PI /180.0)               &
                    + div_x(c2) * sin(vof % q(c2) * PI /180.0)

        div_y(c1) = grid % dy(s) / dotprod * grid % s(s)         &
                    * cos(vof % q(c2) * PI /180.0)               &
                    + div_y(c2) * sin(vof % q(c2) * PI /180.0)

        div_z(c1) = grid % dz(s) / dotprod * grid % s(s)         &
                    * cos(vof % q(c2) * PI /180.0)               &
                    + div_z(c2) * sin(vof % q(c2) * PI /180.0)

        div_x(c2) = div_x(c1)
        div_y(c2) = div_y(c1)
        div_z(c2) = div_z(c1)
      end if

    end if
  end do

  ! Normalize vector at cells
  do c = 1, grid % n_cells
    norm_grad = sqrt(div_x(c) ** 2 + div_y(c) ** 2 + div_z(c) ** 2)
    if (norm_grad >= epsloc) then
      div_x(c) = div_x(c) / norm_grad
      div_y(c) = div_y(c) / norm_grad
      div_z(c) = div_z(c) / norm_grad
    else
      div_x(c) = 0.0
      div_y(c) = 0.0
      div_z(c) = 0.0
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, div_x(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_y(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_z(-nb:nc))

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % curv = 0.0

  ! Derivatives of normals using nodes

  call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
                                              div_x(-nb:nc), grad_x(1:nn),  &
                                           1, div_xx(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_xx(-nb:nc)

  call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
                                              div_y(-nb:nc), grad_y(1:nn),  &
                                           2, div_yy(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_yy(-nb:nc)

  call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
                                              div_z(-nb:nc), grad_z(1:nn),  &
                                           3, div_zz(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_zz(-nb:nc)

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv(-nb:nc))

  call Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,              &
                      div_x(-nb:nc), div_y(-nb:nc), div_z(-nb:nc))

  end subroutine
