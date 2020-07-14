!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Grad_From_Nodes(grid)
!------------------------------------------------------------------------------!
!   Computes nodal weights for gradient calculation at cells using the         !
!   Gram-Schmdt process. This can also be accomplished using T-Flows gradient  !
!   routines                                                                   !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum1 => r_cell_11,  &
                      sum2 => r_cell_12,  &
                      sum3 => r_cell_13,  &
                      r11  => r_cell_14,  &
                      r12  => r_cell_15,  &
                      r13  => r_cell_16,  &
                      r22  => r_cell_17,  &
                      r23  => r_cell_18,  &
                      r33  => r_cell_19
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   ,target :: grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: wx(:,:), wy(:,:), wz(:,:)
  integer                   :: c
  integer                   :: i_nod, n
  !============================================================================!

  ! First take aliases
  wx   => grid % weight_gradx_nodes
  wy   => grid % weight_grady_nodes
  wz   => grid % weight_gradz_nodes

  sum1 = 0.0; sum2 = 0.0; sum3 = 0.0

  do c = 1, grid % n_cells
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)

      sum1(c) = sum1(c) + (grid % xn(n) - grid % xc(c)) ** 2
      sum2(c) = sum2(c) + (grid % xn(n) - grid % xc(c))        &
                        * (grid % yn(n) - grid % yc(c))
      sum3(c) = sum3(c) + (grid % xn(n) - grid % xc(c))        &
                        * (grid % zn(n) - grid % zc(c))
    end do

  end do

  r11 = sqrt(sum1)
  do c = 1, grid % n_cells
    r12(c) = 1.0 / r11(c) * sum2(c)
    r13(c) = 1.0 / r11(c) * sum3(c)
  end do

  sum1 = 0.0; sum2 = 0.0

  do c = 1, grid % n_cells
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)

      sum1(c) = sum1(c) + ( (grid % yn(n) - grid % yc(c))                      &
                        - r12(c) / r11(c) * (grid % xn(n) - grid % xc(c)) ) ** 2
      sum2(c) = sum2(c) + (grid % zn(n) - grid % zc(c))                        &
                        * ( (grid % yn(n) - grid % yc(c))                      &
                        - r12(c) / r11(c) * (grid % xn(n) - grid % xc(c)) )
    end do

  end do


  r22 = sqrt(sum1)

  do c = 1, grid % n_cells
    r23(c) = 1.0 / r22(c) * sum2(c)
  end do

  sum1 = 0.0

  do c = 1, grid % n_cells
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)

      sum1(c) = sum1(c) + ( (grid % zn(n) - grid % zc(c))                      &
                        + (grid % xn(n) - grid % xc(c)) / r11(c)               &
                        * (r23(c) * r12(c) / r22(c) - r13(c))                  &
                        - r23(c)/r22(c) * (grid % yn(n) - grid % yc(c)) ) ** 2
    end do

  end do

  r33 = sqrt(sum1)

  do c = 1, grid % n_cells
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)

      wz(i_nod,c) = 1.0 / r33(c) ** 2 * ( (grid % zn(n) - grid % zc(c))        &
                                        + (grid % xn(n) - grid % xc(c))        &
                                        / r11(c)                               &
                                        * (r23(c) * r12(c) / r22(c) - r13(c))  &
                                        - r23(c)/r22(c)                        &
                                        * (grid % yn(n) - grid % yc(c)) )

      wy(i_nod,c) = 1.0 / r22(c) ** 2 * ( (grid % yn(n) - grid % yc(c))     &
                                        - r12(c) / r11(c)                   &
                                        * (grid % xn(n) - grid % xc(c)) )   &
                                        - r23(c)/r22(c) * wz(i_nod,c)

      wx(i_nod,c) = 1.0 / r11(c) * ( (grid % xn(n) - grid % xc(c)) / r11(c)    &
                                 - r12(c) * wy(i_nod,c)                        &
                                 - r13(c) * wz(i_nod,c) )
    end do

  end do

  end subroutine
