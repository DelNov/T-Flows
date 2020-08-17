!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solve_System(mult, sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type),     target :: sol
  type(Multiphase_Type), target :: mult
  real, contiguous,      target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: vof
  type(Matrix_Type), pointer :: a
  character(SL)              :: solver
  integer                    :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid
  vof => mult % vof
  a   => sol % a

  !! Sum of neighbours
  !neigh(-nb:nc) = 0.0
  !do s = grid % n_bnd_faces + 1, grid % n_faces
  !  c1 = grid % faces_c(1,s)
  !  c2 = grid % faces_c(2,s)
  !  neigh(c1) = neigh(c1) - a % val(a % pos(1,s)) * vof % o(c2)
  !  neigh(c2) = neigh(c2) - a % val(a % pos(2,s)) * vof % o(c1)
  !end do
  !call Grid_Mod_Exchange_Real(grid, neigh(-nb:nc))

  !! Solve velocity explicitely (no under relaxation!!)
  !do c = 1, grid % n_cells
  !  vof % n(c) = (neigh(c) + b(c)) / a % val(a % dia(c))
  !end do

  ! Get solver
  call Control_Mod_Solver_For_Multiphase(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Multiphase')
  call Cg(sol,            &
          vof % n,        &
          b,              &
          vof % precond,  &
          vof % mniter,   &      ! max number of iterations
          vof % eniter,   &      ! executed number of iterations
          vof % tol,      &
          vof % res)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Multiphase')

  !PRINT *, 'VOF ITER ', vof % eniter
!  call Multiphase_Mod_Vof_Scale_Residuals(mult, sol, vof)

  if (.not. mult % phase_change) then
    call Info_Mod_Iter_Fill_At(1, 6, vof % name, vof % eniter, vof % res)
  end if

  end subroutine
