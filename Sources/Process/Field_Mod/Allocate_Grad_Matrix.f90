!==============================================================================!
  subroutine Field_Mod_Allocate_Grad_Matrix(flow)
!------------------------------------------------------------------------------!
!   Allocates memory for gradient matrices                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
!==============================================================================!

  allocate(flow % grad_c2c(6, flow % pnt_grid % n_cells))
  flow % grad_c2c(:,:) = 0.0

  allocate(flow % grad_n2c(6, flow % pnt_grid % n_cells))
  flow % grad_n2c(:,:) = 0.0

  allocate(flow % grad_c2n(6, flow % pnt_grid % n_nodes))
  flow % grad_c2n(:,:) = 0.0

  end subroutine
