!==============================================================================!
  subroutine Var_Mod_Allocate_Solution(name_phi, phi, grid)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)        :: name_phi
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
!==============================================================================!
  
  ! Store variable name
  phi % name = name_phi

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-grid % n_bnd_cells: grid % n_cells));  phi % n  = 0.
  allocate (phi % o (-grid % n_bnd_cells: grid % n_cells));  phi % o  = 0.
  allocate (phi % oo(-grid % n_bnd_cells: grid % n_cells));  phi % oo = 0.

  ! Advection terms
  allocate (phi % a   (grid % n_cells));  phi % a    = 0.
  allocate (phi % a_o (grid % n_cells));  phi % a_o  = 0.
  allocate (phi % a_oo(grid % n_cells));  phi % a_oo = 0.

  ! Diffusion terms
  allocate (phi % d_o (grid % n_cells));  phi % d_o  = 0.
  allocate (phi % d_oo(grid % n_cells));  phi % d_oo = 0.

  ! Cross diffusion terms
  allocate (phi % c  (grid % n_cells));   phi % c    = 0.
  allocate (phi % c_o (grid % n_cells));  phi % c_o  = 0.
  allocate (phi % c_oo(grid % n_cells));  phi % c_oo = 0.

  ! Variable's boundary flux
  allocate (phi % q(-grid % n_bnd_cells: -1));  phi % q  = 0.

  ! Boundary cell type (important for scalars, since they
  ! can have different boundary conditions at the walls)
  allocate (phi % bnd_cell_type(-grid % n_bnd_cells: -1))
  phi % bnd_cell_type = 0

  end subroutine
