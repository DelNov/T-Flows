!==============================================================================!
  subroutine Var_Mod_Create_Variable(phi, Grid, name_phi, name_flux)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi   !! variable object being created
  type(Grid_Type), target :: Grid  !! grid on which it is defined
  character(len=*)        :: name_phi   !! variable's name, connects the
    !! variable to boundary and initial conditions specified in control file
  character(len=*)        :: name_flux  !! name of variable's flux,
    !! connects variable to boundary condition values in control file
!==============================================================================!

  ! Store Grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Variable names must be short, up to 4 (VL) characters long
  Assert(len(name_phi)  .le. VL)
  Assert(len(name_flux) .le. VL)

  ! Store variable name
  phi % name      = name_phi
  phi % flux_name = name_flux

  ! Values in the new (n), old (o) and older than old (oo) time step
  allocate(phi % n (-Grid % n_bnd_cells:Grid % n_cells));  phi % n  = 0.0
  allocate(phi % o (-Grid % n_bnd_cells:Grid % n_cells));  phi % o  = 0.0
  if(phi % td_scheme .eq. PARABOLIC) then
    allocate(phi % oo(-Grid % n_bnd_cells:Grid % n_cells));  phi % oo = 0.0
  end if

  ! Variable's boundary value
  allocate (phi % b(-Grid % n_bnd_cells:-1));  phi % b = 0.0

  ! Variable's boundary flux
  ! (It has to go through inside cells to be able to call Exchange on it.)
  allocate (phi % q(-Grid % n_bnd_cells:Grid % n_cells));  phi % q = 0.0

  ! Boundary cell type (important for scalars, since they
  ! can have different boundary conditions at the walls)
  ! It expands over all faces, in the case we decide to store ...
  ! ... store boundary conditions in the faces one day, and ...
  ! ... thus get rid of the "if(c2 < 0) then" checks
  allocate (phi % bnd_cond_type(-Grid % n_bnd_cells : Grid % n_faces))
  phi % bnd_cond_type = 0

  end subroutine
