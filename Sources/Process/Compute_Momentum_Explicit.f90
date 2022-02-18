!==============================================================================!
  subroutine Compute_Momentum_Explicit(Flow, ui, Nat)
!------------------------------------------------------------------------------!
!   Explicit computation of momentum equations, used in PISO algorithm,        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Numerics_Mod
  use Work_Mod,     only: neigh => r_cell_01
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                    (allocates Work_Mod)           !
!     |                                                                        !
!     +----> Compute_Momentum                   (doesn't use Work_Mod)         !
!             |                                                                !
!             +----> Compute_Momentum_Explicit  (safe to use r_cell_01)        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Var_Type)            :: ui        ! velocity component
  type(Native_Type), target :: Nat
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Matrix_Type), pointer :: M
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  M    => Nat % M
  b    => Nat % b % val

  ! PISO corrections are executed here
  if (Flow % p_m_coupling == PISO .and.  &
      Flow % piso_status .eqv. .true.) then

    ! Sum of neighbours
    neigh = 0.0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 > 0) then
        neigh(c1) = neigh(c1) - M % val(M % pos(1,s)) * ui % n(c2)
        neigh(c2) = neigh(c2) - M % val(M % pos(2,s)) * ui % n(c1)
      end if
    end do
    call Grid % Exchange_Cells_Real(neigh)

    ! Solve velocity explicitely (no under relaxation!!)
    do c = 1, Grid % n_cells
      ui % n(c) = (neigh(c) + b(c)) / M % val(M % dia(c))
    end do

    call Flow % Grad_Variable(ui)

  end if

  end subroutine
