!==============================================================================!
  subroutine Create_Field(Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
  type(Grid_Type),   target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: nb, nc, ns
!==============================================================================!

  ! Store the pointer to a Grid
  Flow % pnt_grid => Grid

  ! Take some aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
  ns = Grid % n_faces

  !--------------------------------------------------------------!
  !   Create native solvers (matrices A, M, right hand side b,   !
  !   helping vectors for CG method such as p, q, r and d_inv)   !
  !--------------------------------------------------------------!
  print '(a)', ' # Creating field''s native solver'
  call Flow % Nat % Create_Native(Grid)

  !----------------------------------!
  !   Memory for gradient matrices   !
  !----------------------------------!
  allocate(Flow % grad_c2c(6, nc));  Flow % grad_c2c(:,:) = 0.0

  !------------------------!
  !   Memory for unknows   !
  !------------------------!

  ! Momentum equations
  call Var_Mod_Create_Variable(Flow % u, Grid)
  call Var_Mod_Create_Variable(Flow % v, Grid)
  call Var_Mod_Create_Variable(Flow % w, Grid)

  ! Pressure
  call Var_Mod_Create_Variable(Flow % pp, Grid)  ! pressure correction
  call Var_Mod_Create_Variable(Flow % p,  Grid)  ! pressure

  ! Volume fluxes through faces
  allocate(Flow % v_flux(ns))

  end subroutine
