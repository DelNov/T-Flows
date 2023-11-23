!==============================================================================!
  subroutine Var_Mod_Create_Solution(phi, A, name_phi, name_flux,  &
                                     reuse_pet)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the Grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)             :: phi
  type(Matrix_Type),  target :: A
  character(len=*)           :: name_phi
  character(len=*)           :: name_flux
  integer,          optional :: reuse_pet
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
