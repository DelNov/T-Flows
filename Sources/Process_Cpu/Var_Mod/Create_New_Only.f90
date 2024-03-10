!==============================================================================!
  subroutine Var_Mod_Create_New_Only(phi, Grid, name_phi)
!------------------------------------------------------------------------------!
!>  This subroutine allocates a simplified unknown variable, which holds only
!>  its current value (n) and graduents (x, y, z).  It's apt for variables
!>  like pressure or smoothed variant of VOF, which are calculated from other
!>  variables' values, either by accumulation or convolution (smoothing).
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Associates the variable (phi) with its grid (Grid) and sets its name     !
!     (name_phi). It doesn't link to a matrix as the variable isn't solved     !
!     through a PDE.                                                           !
!   * Simplification: Unlike Var_Mod_Create_Solution, this subroutine only     !
!     allocates and initializes the current value and the variable's gradient  !
!     and boundary values. It omits aspects like old values and linear solver  !
!     parameters, as they are not needed for this type of variable.            !
!   * Usage: Ideal for variables that require less complexity in terms of      !
!     storage and computation, focusing on current values, gradients, and      !
!     boundary conditions.                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi        !! variable object being created
  type(Grid_Type), target :: Grid       !! grid on which it is defined
  character(len=*)        :: name_phi   !! variable's name, connects the
    !! variable to boundary and initial conditions specified in control file
!==============================================================================!

  ! Store Grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n(-Grid % n_bnd_cells:Grid % n_cells));  phi % n = 0.0

  ! Gradients
  allocate (phi % x(-Grid % n_bnd_cells:Grid % n_cells));  phi % x = 0.0
  allocate (phi % y(-Grid % n_bnd_cells:Grid % n_cells));  phi % y = 0.0
  allocate (phi % z(-Grid % n_bnd_cells:Grid % n_cells));  phi % z = 0.0

  ! Variable's boundary value
  allocate (phi % b(-Grid % n_bnd_cells:Grid % n_cells));  phi % b = 0.

  end subroutine
