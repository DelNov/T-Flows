!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Nodal_Grad(grid)
!------------------------------------------------------------------------------!
!   Computes cell weights for gradient calculation at nodes using the          !
!   Gram-Schmdt process.                                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: sum1   => r_node_01,  &
                      sum2   => r_node_02,  &
                      sum3   => r_node_03,  &
                      r11    => r_node_04,  &
                      r12    => r_node_05,  &
                      r13    => r_node_06,  &
                      r22    => r_node_07,  &
                      r23    => r_node_08,  &
                      r33    => r_node_09
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   ,target :: grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: wx(:,:), wy(:,:), wz(:,:)
  integer                   :: c
  integer                   :: i_nod, i_cell, n
  real                      :: epsloc
  !============================================================================!

  epsloc = epsilon(epsloc)

  ! First take aliases
  allocate(grid % weight_gradx_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))
  allocate(grid % weight_grady_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))
  allocate(grid % weight_gradz_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))

  wx   => grid % weight_gradx_cells
  wy   => grid % weight_grady_cells
  wz   => grid % weight_gradz_cells

  sum1 = 0.0; sum2 = 0.0; sum3 = 0.0

  do n = 1, grid % n_nodes
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)

      sum1(n) = sum1(n) + (grid % xc(c) - grid % xn(n)) ** 2
      sum2(n) = sum2(n) + (grid % xc(c) - grid % xn(n))        &
                        * (grid % yc(c) - grid % yn(n))
      sum3(n) = sum3(n) + (grid % xc(c) - grid % xn(n))        &
                        * (grid % zc(c) - grid % zn(n))
    end do

  end do

  r11 = sqrt(sum1)
  do n = 1, grid % n_nodes
    r12(n) = 1.0 / (r11(n) + epsloc) * sum2(n)
    r13(n) = 1.0 / (r11(n) + epsloc) * sum3(n)
  end do

  sum1 = 0.0; sum2 = 0.0

  do n = 1, grid % n_nodes
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)

      sum1(n) = sum1(n) + ( (grid % yc(c) - grid % yn(n))                      &
                        - r12(n) / r11(n) * (grid % xc(c) - grid % xn(n)) ) ** 2
      sum2(n) = sum2(n) + (grid % zc(c) - grid % zn(n))                        &
                        * ( (grid % yc(c) - grid % yn(n))                      &
                        - r12(n) / r11(n) * (grid % xc(c) - grid % xn(n)) )
    end do

  end do


  r22 = sqrt(sum1)

  do n = 1, grid % n_nodes
    r23(n) = 1.0 / (r22(n) + epsloc) * sum2(n)
  end do

  sum1 = 0.0

  do n = 1, grid % n_nodes
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)

      sum1(n) = sum1(n) + ( (grid % zc(c) - grid % zn(n))                      &
                        + (grid % xc(c) - grid % xn(n)) / r11(n)               &
                        * (r23(n) * r12(n) / (r22(n) + epsloc) - r13(n))       &
                        - r23(n)/(r22(n) + epsloc)                             &
                        * (grid % yc(c) - grid % yn(n)) ) ** 2
    end do

  end do

  r33 = sqrt(sum1)

  do n = 1, grid % n_nodes
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)

      wz(i_cell,n) = 1.0 / (r33(n) ** 2 + epsloc)    &
                                       * ( (grid % zc(c) - grid % zn(n))       &
                                         + (grid % xc(c) - grid % xn(n))       &
                                         / (r11(n) + epsloc)                   &
                                         * (r23(n) * r12(n) / (r22(n) + epsloc)&
                                         - r13(n))                             &
                                         - r23(n)/(r22(n) + epsloc)            &
                                         * (grid % yc(c) - grid % yn(n)) )

      wy(i_cell,n) = 1.0 / (r22(n) ** 2 + epsloc)    &
                                        * ( (grid % yc(c) - grid % yn(n))     &
                                         - r12(n) / r11(n)                   &
                                         * (grid % xc(c) - grid % xn(n)) )   &
                                         - r23(n)/(r22(n) + epsloc)          &
                                         * wz(i_cell,n)

      wx(i_cell,n) = 1.0 / r11(n) * ( (grid % xc(c) - grid % xn(n)) / r11(n)   &
                                  - r12(n) * wy(i_cell,n)                      &
                                  - r13(n) * wz(i_cell,n) )
    end do

  end do

  end subroutine
