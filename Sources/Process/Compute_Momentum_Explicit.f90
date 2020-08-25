!==============================================================================!
  subroutine Compute_Momentum_Explicit(flow, i, sol)
!------------------------------------------------------------------------------!
!   Explicit computation of momentum equations, used in PISO algorithm,        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: neigh => r_cell_12,   &
                      res   => r_cell_13
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  integer                       :: i           ! component
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: ui
  type(Matrix_Type),pointer :: a
  real, contiguous, pointer :: b(:)
  integer                   :: s, c, c1, c2, nt, ni
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  if(i .eq. 1) ui => flow % u
  if(i .eq. 2) ui => flow % v
  if(i .eq. 3) ui => flow % w

  call Solver_Mod_Alias_System(sol, a, b)

  ! PISO corrections are executed here
  if (flow % p_m_coupling == PISO .and.  &
      flow % piso_status .eqv. .true.) then

    ! Sum of neighbours
    neigh = 0.0
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      neigh(c1) = neigh(c1) - a % val(a % pos(1,s)) * ui % n(c2)
      neigh(c2) = neigh(c2) - a % val(a % pos(2,s)) * ui % n(c1)
    end do
    call Grid_Mod_Exchange_Cells_Real(grid, neigh)

    ! Solve velocity explicitely (no under relaxation!!)
    do c = 1, grid % n_cells
      ui % n(c) = (neigh(c) + b(c)) / a % val(a % dia(c))
    end do

    call Field_Mod_Grad_Variable(flow, ui)

    if (flow % i_corr == flow % n_piso_corrections) then
      res(:) = 0.0
      nt = a % pnt_grid % n_cells
      ni = a % pnt_grid % n_cells - a % pnt_grid % comm % n_buff_cells
      call Residual_Vector(ni, res(1:nt), b(1:nt), a, ui % n(1:nt))
      ui % res_scal = sum(abs(res(1:ni)))
      call Comm_Mod_Global_Sum_Real(ui % res_scal)
    end if

  end if

  end subroutine
