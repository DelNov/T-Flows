!==============================================================================!
  subroutine Grad_Modeled_Flow(Swarm)
!------------------------------------------------------------------------------!
!>  Calculates and stores the gradients of modeled flow parameters, namely
!>  the v2 component from the k-eps-zeta-f model. It computes gradients of
!>  turbulent kinetic energy and dissipation rate, althogh I am not sure it
!>  uses them. It computes v2_mod from zeta and turbulent kinetic energy and
!>  then computes its gradients in x, y and z direction.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Turbulence data access: Gathers turbulence data from the flow field,     !
!     grid, and turbulence model for gradient calculations.                    !
!   * Gradient calculation: Computes spatial gradients of turbulent kinetic    !
!     energy and dissipation rate, although I am not sure it is needed.        !!
!   * Modeled turbulence parameter: Forms an array of a modeled turbulence     !
!     parameter (v^2) based on turbulence data.                                !
!   * Gradient storage: Stores spatial gradients of the modeled turbulence     !
!     variable v2 in three dimensions, important for particle dynamics.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Turb_Type),  pointer :: Turb
  integer                   :: c
!==============================================================================!

  ! Take aliases for Flow
  Flow => Swarm % pnt_flow
  Grid => Swarm % pnt_grid
  Turb => Swarm % pnt_turb

  ! Gradients of turbulent quantities
  call Flow % Grad_Variable(Turb % kin)
  call Flow % Grad_Variable(Turb % eps)

  ! Array for v^2
  do c = Cells_In_Domain_And_Buffers()
    Swarm % v2_mod(c) = sqrt(Turb % zeta % n(c) * Turb % kin % n(c))
  end do

  ! Storing the gradients of v2_mod in Work_Mod arrays
  ! (dv^2/dx, dv^2/dy, dv^2/dz)
  call Flow % Grad_Component(Grid, Swarm % v2_mod, 1, Swarm % v2_mod_x)
  call Flow % Grad_Component(Grid, Swarm % v2_mod, 2, Swarm % v2_mod_y)
  call Flow % Grad_Component(Grid, Swarm % v2_mod, 3, Swarm % v2_mod_z)

  end subroutine
