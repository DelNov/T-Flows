!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Csf(mult)
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
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: smooth
  integer                   :: c, c1, c2, s, nb, nc
  real                      :: v1(3), v2(3), v3(3), v4(3)
  real                      :: norm_grad, epsloc, dotprod
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
    norm_grad = sqrt(  smooth % x(c) ** 2  &
                     + smooth % y(c) ** 2  &
                     + smooth % z(c) ** 2)
    if(norm_grad >= epsloc) then
      mult % nx(c) = smooth % x(c) / norm_grad
      mult % ny(c) = smooth % y(c) / norm_grad
      mult % nz(c) = smooth % z(c) / norm_grad
    else
      mult % nx(c) = 0.0
      mult % ny(c) = 0.0
      mult % nz(c) = 0.0
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, mult % nx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, mult % ny(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, mult % nz(-nb:nc))

  ! Tangent vector to walls/symmetries

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.   &
        Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then

      norm_grad = norm2((/mult % nx(c1), mult % ny(c1), mult % nz(c1)/))

      if(norm_grad > epsloc) then
        v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
        v2 = (/mult % nx(c1), mult % ny(c1), mult % nz(c1)/)
        v3 = Math_Mod_Cross_Product(v1, v2)
        v4 = Math_Mod_Cross_Product(v3, v1)
        ! projection on v4
        norm_grad = norm2(v4)
        if(norm_grad > epsloc) then
          mult % nx(c2) = v4(1) / norm_grad
          mult % ny(c2) = v4(2) / norm_grad
          mult % nz(c2) = v4(3) / norm_grad
        end if
      end if
    else
      mult % nx(c2) = mult % nx(c1)
      mult % ny(c2) = mult % ny(c1)
      mult % nz(c2) = mult % nz(c1)
    end if
  end do

  ! Correct for contact angle at walls

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then

      ! Accumulate values of faces
      norm_grad = norm2((/mult % nx(c1),mult % ny(c1),mult % nz(c1)/))
      if(norm_grad > epsloc) then
        dotprod = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/),   &
                              (/grid % sx(s), grid % sy(s), grid % sz(s)/))

        mult % nx(c1) = grid % dx(s) / dotprod * grid % s(s)                  &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + mult % nx(c2) * sin(vof % q(c2) * PI /180.0)

        mult % ny(c1) = grid % dy(s) / dotprod * grid % s(s)                  &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + mult % ny(c2) * sin(vof % q(c2) * PI /180.0)

        mult % nz(c1) = grid % dz(s) / dotprod * grid % s(s)                  &
                    * cos(vof % q(c2) * PI /180.0)                            &
                    + mult % nz(c2) * sin(vof % q(c2) * PI /180.0)

        mult % nx(c2) = mult % nx(c1)
        mult % ny(c2) = mult % ny(c1)
        mult % nz(c2) = mult % nz(c1)
      end if

    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, mult % nx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, mult % ny(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, mult % nz(-nb:nc))

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  ! Find divergence of normals
  call Field_Mod_Grad_Component(flow, mult % nx(-nb:nc), 1, div_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, mult % ny(-nb:nc), 2, div_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, mult % nz(-nb:nc), 3, div_z(-nb:nc))

  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_x(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_y(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_z(-nb:nc)

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

! call Grid_Mod_Save_Debug_Vtu(grid, 'curv_sharp', scalar_cell = mult % curv,  &
!                                                  scalar_name = 'curv_sharp')

  call Multiphase_Mod_Vof_Smooth_Curvature(mult)

! call Grid_Mod_Save_Debug_Vtu(grid, 'curv_smooth', scalar_cell = mult % curv,  &
!                                                   scalar_name = 'curv_smooth')

  do c = 1, grid % n_cells
    if(smooth % n(c) < 0.01 .or. smooth % n(c) > 0.99) mult % curv(c) = 0.0
  end do

  end subroutine
