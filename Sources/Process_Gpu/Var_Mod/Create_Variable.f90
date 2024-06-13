!==============================================================================!
  subroutine Var_Mod_Create_Variable(phi, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi   !! variable object being created
  type(Grid_Type), target :: Grid  !! grid on which it is defined
!==============================================================================!

  ! Store Grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Values in the new (n) and the old (o) time step
  allocate (phi % n(-Grid % n_bnd_cells:Grid % n_cells));  phi % n = 0.0
  allocate (phi % o(-Grid % n_bnd_cells:Grid % n_cells));  phi % o = 0.0

  ! Variable's boundary value
  allocate (phi % b(-Grid % n_bnd_cells:-1));  phi % b = 0.

  ! Variable's boundary flux
  allocate (phi % q(-Grid % n_bnd_cells:-1));  phi % q = 0.

  ! Boundary cell type (important for scalars, since they
  ! can have different boundary conditions at the walls)
  ! It expands over all faces, in the case we decide to store ...
  ! ... store boundary conditions in the faces one day, and ...
  ! ... thus get rid of the "if(c2 < 0) then" checks
  allocate (phi % bnd_cond_type(-Grid % n_bnd_cells : Grid % n_faces))
  phi % bnd_cond_type = 0

  ! Gradient components
  allocate (phi % x(-Grid % n_bnd_cells:Grid % n_cells));  phi % x = 0.0
  allocate (phi % y(-Grid % n_bnd_cells:Grid % n_cells));  phi % y = 0.0
  allocate (phi % z(-Grid % n_bnd_cells:Grid % n_cells));  phi % z = 0.0

  end subroutine
