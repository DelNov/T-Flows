!==============================================================================!
  subroutine Calculate_Vorticity(flow)
!------------------------------------------------------------------------------!
!  Computes the magnitude of the vorticity                                     !
!------------------------------------------------------------------------------!
!  vort = sqrt( 2 * Sij * Sij )                                                !
!  Sij = 1/2 ( dUi/dXj - dUj/dXi )                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  call Comm_Mod_Exchange_Real(grid, u % n)
  call Comm_Mod_Exchange_Real(grid, v % n)
  call Comm_Mod_Exchange_Real(grid, w % n)

  call Grad_Mod_For_Phi(grid, u % n, 1, u % x, .true.)  ! du/dx
  call Grad_Mod_For_Phi(grid, u % n, 2, u % y, .true.)  ! du/dy
  call Grad_Mod_For_Phi(grid, u % n, 3, u % z, .true.)  ! du/dz

  call Grad_Mod_For_Phi(grid, v % n, 1, v % x, .true.)  ! dv/dx
  call Grad_Mod_For_Phi(grid, v % n, 2, v % y, .true.)  ! dv/dy
  call Grad_Mod_For_Phi(grid, v % n, 3, v % z, .true.)  ! dv/dz

  call Grad_Mod_For_Phi(grid, w % n, 1, w % x, .true.)  ! dw/dx
  call Grad_Mod_For_Phi(grid, w % n, 2, w % y, .true.)  ! dw/dy
  call Grad_Mod_For_Phi(grid, w % n, 3, w % z, .true.)  ! dw/dz

  do c = 1, grid % n_cells
    vort(c) = 2.0 * (0.5 * (w % y(c) - v % z(c)))**2  &
            + 2.0 * (0.5 * (w % x(c) - u % z(c)))**2  &
            + 2.0 * (0.5 * (v % x(c) - u % y(c)))**2

    vort(c) = sqrt(abs(2.0 * vort(c)))
  end do

  end subroutine
