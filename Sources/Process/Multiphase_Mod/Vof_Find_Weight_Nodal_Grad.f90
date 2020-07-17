!==============================================================================!
  subroutine Multiphase_Mod_Vof_Find_Weight_Nodal_Grad(grid)
!------------------------------------------------------------------------------!
!   Computes cell weights for gradient calculation at nodes using the          !
!   Gram-Schmdt process.                                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: wx(:,:), wy(:,:), wz(:,:)
  real, allocatable         :: sum1(:)
  real, allocatable         :: sum2(:)
  real, allocatable         :: sum3(:)
  real, allocatable         :: r11 (:)
  real, allocatable         :: r12 (:)
  real, allocatable         :: r13 (:)
  real, allocatable         :: r22 (:)
  real, allocatable         :: r23 (:)
  real, allocatable         :: r33 (:)
  integer                   :: c
  integer                   :: i_nod, i_cell, n
  real                      :: epsloc
!==============================================================================!

  ! Allocate local arrays
  allocate(sum1(1 : grid % n_nodes))
  allocate(sum2(1 : grid % n_nodes))
  allocate(sum3(1 : grid % n_nodes))
  allocate(r11 (1 : grid % n_nodes))
  allocate(r12 (1 : grid % n_nodes))
  allocate(r13 (1 : grid % n_nodes))
  allocate(r22 (1 : grid % n_nodes))
  allocate(r23 (1 : grid % n_nodes))
  allocate(r33 (1 : grid % n_nodes))

  epsloc = epsilon(epsloc)

  allocate(grid % weight_gradx_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))
  allocate(grid % weight_grady_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))
  allocate(grid % weight_gradz_cells(size(grid % nodes_c,1),   &
                                     size(grid % nodes_c,2)))

  ! Take aliases
  wx => grid % weight_gradx_cells
  wy => grid % weight_grady_cells
  wz => grid % weight_gradz_cells

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

  ! Deallocate local arrays
  deallocate(sum1)
  deallocate(sum2)
  deallocate(sum3)
  deallocate(r11 )
  deallocate(r12 )
  deallocate(r13 )
  deallocate(r22 )
  deallocate(r23 )
  deallocate(r33 )

  end subroutine
