!==============================================================================!
  subroutine Numerics_Mod_Inertial_Term(phi, coef, sol, dt)
!------------------------------------------------------------------------------!
!   Purpose: Dicretize inertial term in conservation equations.                !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Var_Type)            :: phi
  real                      :: coef(-phi % pnt_grid % n_bnd_cells:  &
                                     phi % pnt_grid % n_cells)
  type(Solver_Type), target :: sol
  real                      :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  real                       :: a0
  integer                    :: c
!==============================================================================!

  ! Take alias to grid
  grid => phi % pnt_grid
  a    => sol % a
  b    => sol % b % val

  ! Two time levels; Linear interpolation
  if(phi % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = coef(c) * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(phi % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = coef(c) * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  end subroutine
