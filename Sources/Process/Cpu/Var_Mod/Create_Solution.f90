!==============================================================================!
  subroutine Var_Mod_Create_Solution(phi, A, name_phi, name_flux,  &
                                     reuse_pet)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to allocate and initialize a Var_Type variable
!>  whose values are obtained from a solution of partial diffetential equations.
!>  Examples of such variable include velocity components, temperature,
!>  variables used in turbulence modeling (k, eps, f22, ...), VOF, pressure
!>  correction and even wall distance or potential for velocity initialization.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * It links the variable (phi) to its corresponding grid (Grid) and matrix  !
!     (A), indicating the structure and relationships necessary for the solver.!
!   * PETSc Integration: If PETSc is enabled (conditional compilation), the    !
!     subroutine either creates a new PETSc instance for the variable or       !
!     reuses an existing one, depending on the reuse_pet parameter. This       !
!     unique PETSc instance per variable approach is consistent with the       !
!     overall design philosophy of Var_Mod.                                    !
!   * Initialization: Allocates memory for the variable's current, old, and    !
!     older-than-old values (n, o, oo), as well as for boundary values,        !
!     fluxes, and gradient components. It also initializes these arrays to     !
!     zero and sets up boundary condition types.                               !
!   * Flexibility: The subroutine offers flexibility in managing the linear    !
!     solver settings for each variable, aligning with the module's approach   !
!     to handle each variable's solver requirements individually.              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)             :: phi        !! variable object being created
  type(Matrix_Type),  target :: A          !! system matrix to
  character(len=*)           :: name_phi   !! variable's name, connects the
    !! variable to boundary and initial conditions specified in control file
  character(len=*)           :: name_flux  !! name of variable's flux,
    !! connects variable to boundary condition values in control file
  integer,          optional :: reuse_pet
    !! if true, existing PETSc instance, useful for velocity components
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Store Grid for which the variable is defined
  Grid             => A % pnt_grid
  phi % pnt_grid   => Grid
  phi % pnt_matrix => A

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

# if T_FLOWS_PETSC == 1
  ! Variable creates its own new PETSc type
  if(.not. present(reuse_pet)) then
    Work_Pet % n_members = Work_Pet % n_members + 1
    call Work_Pet % Member(Work_Pet % n_members) % Create_Petsc(  &
      phi % pnt_matrix,                                           &
      phi % name,                                                 &
      phi % pet_rank                                              &
    )
    Assert(phi % pet_rank .eq. Work_Pet % n_members)
    phi % Pet => Work_Pet % Member(Work_Pet % n_members)

  ! Variable re-uses the existing
  else
    Assert(reuse_pet > 0)
    Assert(reuse_pet <= Work_Pet % n_members)
    phi % pet_rank = reuse_pet
    phi % Pet => Work_Pet % Member(reuse_pet)
  end if
# else
  Unused(reuse_pet)
# endif

  end subroutine
