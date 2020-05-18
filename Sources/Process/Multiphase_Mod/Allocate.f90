!==============================================================================!
  subroutine Multiphase_Mod_Allocate(mult, flow)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Field_Type),      target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: nb, nc, nf
!==============================================================================!

  ! Store pointers
  mult % pnt_flow => flow
  mult % pnt_grid => flow % pnt_grid

  ! Take aliases
  grid => flow % pnt_grid
  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nf = grid % n_faces

  call Var_Mod_Allocate_Solution(mult % vof, grid, 'VOF', '')

  if (mult % d_func) then
    call Var_Mod_Allocate_Solution(mult % dist_func, grid, 'D_FUNC', '')
  end if

  allocate(mult % vof_f(nf));  mult % vof_f(1:nf) = 0.0

  allocate(mult % curv(-nb:nc));  mult % curv(-nb:nc) = 0.0

  ! Physical properties for all (two) phases
  allocate(mult % phase_dens(2))
  allocate(mult % phase_visc(2))

  end subroutine
