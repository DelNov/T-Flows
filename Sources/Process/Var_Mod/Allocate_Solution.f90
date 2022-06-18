!==============================================================================!
  subroutine Var_Mod_Allocate_Solution(phi, Grid, name_phi, name_flux)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the Grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi
  type(Grid_Type), target :: Grid
  character(len=*)        :: name_phi
  character(len=*)        :: name_flux
!==============================================================================!

  ! Store Grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Store variable name
  phi % name      = name_phi
  phi % flux_name = name_flux

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-Grid % n_bnd_cells:Grid % n_cells));  phi % n  = 0.
  allocate (phi % o (-Grid % n_bnd_cells:Grid % n_cells));  phi % o  = 0.
  allocate (phi % oo(-Grid % n_bnd_cells:Grid % n_cells));  phi % oo = 0.

  ! Variable's boundary value
  allocate (phi % b(-Grid % n_bnd_cells:Grid % n_cells));  phi % b = 0.

  ! Variable's boundary flux
  allocate (phi % q(-Grid % n_bnd_cells:Grid % n_cells));  phi % q = 0.

  ! Boundary cell type (important for scalars, since they
  ! can have different boundary conditions at the walls)
  ! It expands over all faces, in the case we decide to store ...
  ! ... store boundary conditions in the faces one day, and ...
  ! ... thus get rid of the "if(c2 < 0) then" checks
  allocate (phi % bnd_cond_type(-Grid % n_bnd_cells : Grid % n_faces))
  phi % bnd_cond_type = 0

  ! Gradients
  allocate (phi % x(-Grid % n_bnd_cells:Grid % n_cells));  phi % x = 0.0
  allocate (phi % y(-Grid % n_bnd_cells:Grid % n_cells));  phi % y = 0.0
  allocate (phi % z(-Grid % n_bnd_cells:Grid % n_cells));  phi % z = 0.0

  end subroutine
