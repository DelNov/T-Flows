!==============================================================================!
  subroutine Multiphase_Mod_Vof_Heaviside_Function(mult)
!------------------------------------------------------------------------------!
!   Computes the Heaviside function, necessary for Brackbill's CSF approach    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: heaviside => r_cell_13
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: dist_func
  integer                   :: c, s, c1, c2, nb, nc
  real                      :: eps_grid
!==============================================================================!

  flow      => mult % pnt_flow
  grid      => mult % pnt_grid
  vof       => mult % vof
  dist_func => mult % dist_func

  nb = grid % n_bnd_cells
  nc = grid % n_cells


  do c = 1, grid % n_cells

    eps_grid = 1.5 * grid % vol(c) ** ONE_THIRD

    if (dist_func % n(c) > eps_grid) then
      heaviside(c) = 1.0
    else if (dist_func % n(c) < -eps_grid) then
      heaviside(c) = 0.0
    else
      heaviside(c) = 0.5 * ( 1.0 + dist_func % n(c) / eps_grid       &
                            + 1.0 / PI * sin(PI * dist_func % n(c)    &
                            / eps_grid))
    end if

  end do

  ! Values at boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (c2 < 0) then
      heaviside(c2) = heaviside(c1)
    end if
  end do

  ! Find gradients
  call Field_Mod_Grad_Component(flow, heaviside(-nb:nc), 1, vof % x)
  call Field_Mod_Grad_Component(flow, heaviside(-nb:nc), 2, vof % y)
  call Field_Mod_Grad_Component(flow, heaviside(-nb:nc), 3, vof % z)

  ! Store Heavyside function for the pressure correction step
  dist_func % oo(-nb:nc) = heaviside(-nb:nc)

  end subroutine
