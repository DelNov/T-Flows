!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(grid)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the shear stress.                                !
!------------------------------------------------------------------------------!
!   shear = sqrt( 2 * Sij * Sij )                                              !
!   Sij = 1/2 ( dUi/dXj + dUj/dXi )                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grad_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!
  
  call Comm_Mod_Exchange(grid, u % n)
  call Comm_Mod_Exchange(grid, v % n)
  call Comm_Mod_Exchange(grid, w % n)

  !---------------!
  !   SGS terms   !
  !---------------!
  call Grad_Mod_For_Phi(grid, u % n, 1, u % x, .true.)  ! du/dx
  call Grad_Mod_For_Phi(grid, u % n, 2, u % y, .true.)  ! du/dy
  call Grad_Mod_For_Phi(grid, u % n, 3, u % z, .true.)  ! du/dz
  call Grad_Mod_For_Phi(grid, v % n, 1, v % x, .true.)  ! dv/dx
  call Grad_Mod_For_Phi(grid, v % n, 2, v % y, .true.)  ! dv/dy
  call Grad_Mod_For_Phi(grid, v % n, 3, v % z, .true.)  ! dv/dz
  call Grad_Mod_For_Phi(grid, w % n, 1, w % x, .true.)  ! dw/dx
  call Grad_Mod_For_Phi(grid, w % n, 2, w % y, .true.)  ! dw/dy
  call Grad_Mod_For_Phi(grid, w % n, 3, w % z, .true.)  ! dw/dz

  shear(:) =  u % x(:)**2.                    &
            + v % y(:)**2.                    &
            + w % z(:)**2.                    &
            + 0.5 * (v % z(:) + w % y(:))**2. & 
            + 0.5 * (u % z(:) + w % x(:))**2. & 
            + 0.5 * (v % x(:) + u % y(:))**2.

  vort(:) = - (   0.5*(v % z(:) - w % y(:))**2. &
                + 0.5*(u % z(:) - w % x(:))**2. &
                + 0.5*(v % x(:) - u % y(:))**2.)

  shear(:) = sqrt(2.0 * shear(:))
  vort(:)  = sqrt(2.0 * abs(vort(:)))

  end subroutine
