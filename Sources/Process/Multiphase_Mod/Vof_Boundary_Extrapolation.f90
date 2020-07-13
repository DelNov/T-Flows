!==============================================================================!
  subroutine Multiphase_Mod_Vof_Boundary_Extrapolation(grid, mult, var)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Gauss theorem        !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: grad_x => r_cell_21,   &
                      grad_y => r_cell_22,   &
                      grad_z => r_cell_23
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)              :: grid
  type(Multiphase_Type), target :: mult
  real                          :: var(-grid % n_bnd_cells    &
                                      : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  integer                       :: s, ss, c, c1, c2, n
  integer                       :: i_fac, count_neg, cc1, cc2, cc
  real                          :: gx, gy, gz
  real                          :: epsloc, num, var_f, fss, signo
!==============================================================================!

  vof  => mult % vof
  flow => mult % pnt_flow

  epsloc = epsilon(epsloc)

  ! Simply zero gradient at boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
      var(c2) = var(c1)
    end if
  end do

  ! A preliminar gradient
  call Multiphase_Mod_Vof_Grad_Component(flow, var, 1, grad_x)
  call Multiphase_Mod_Vof_Grad_Component(flow, var, 2, grad_y)
  call Multiphase_Mod_Vof_Grad_Component(flow, var, 3, grad_z)

  ! First step: estimate face values for all faces at wall boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then

      count_neg = 0
      ! Loop on faces
      gx = 0.0; gy = 0.0; gz = 0.0
      do i_fac = 1, grid % cells_n_faces(c1)
        ss = grid % cells_f(i_fac, c1)
        if (s .ne. ss) then !average gradients
          cc1 = grid % faces_c(1,ss)
          cc2 = grid % faces_c(2,ss)
          if (c1 == cc1) then
            cc = cc2
          else
            cc = cc1
          end if
          gx = gx + grad_x(cc)
          gy = gy + grad_y(cc)
          gz = gz + grad_z(cc)
        end if
      end do

      gx = gx / real(grid % cells_n_faces(c1))
      gy = gy / real(grid % cells_n_faces(c1))
      gz = gz / real(grid % cells_n_faces(c1))

      var(c2) = var(c1)                                                     &
              + dot_product((/gx, gy, gz/),                                 &
                            (/grid % dx(s), grid % dy(s), grid % dz(s)/))
    end if
  end do

  ! With divergence theorem:
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then

      count_neg = 0
      ! Loop on faces
      gx = 0.0; gy = 0.0; gz = 0.0
      do i_fac = 1, grid % cells_n_faces(c1)
        ss = grid % cells_f(i_fac, c1)
        fss = grid % fw(ss)
        num = 0.0

        if (s .ne. ss) then !average gradients
          cc1 = grid % faces_c(1,ss)
          cc2 = grid % faces_c(2,ss)
          if (c1 == cc1) then
            signo = 1.0
          else
            signo = -1.0
          end if

          var_f = fss * var(cc1) + (1.0 - fss) * var(c2)
          num = num + var_f * signo * dot_product(                            &
                              (/grid % dx(ss), grid % dy(ss), grid % dz(ss)/),   &
                              (/grid % sx(ss), grid % sy(ss), grid % sz(ss)/))
        end if
      end do

      num = num / grid % vol(c1)

      var(c2) = (var(c1) + num) / (1.0 - dot_product(                         &
                               (/grid % dx(s), grid % dy(s), grid % dz(s)/),  &
                               (/grid % sx(s), grid % sy(s), grid % sz(s)/))  &
                                /grid % vol(c1))
    end if
  end do
  end subroutine
