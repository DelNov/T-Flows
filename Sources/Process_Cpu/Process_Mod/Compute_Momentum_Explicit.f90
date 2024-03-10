!==============================================================================!
  subroutine Compute_Momentum_Explicit(Process, Flow, ui, Nat)
!------------------------------------------------------------------------------!
!>  This subroutine explicitly computes the momentum equations, particularly
!>  used within the PISO algorithm framework.  It is integral in updating
!>  velocity fields based on momentum conservation.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and setup: Initiates performance profiling and sets up    !
!     required local variables and pointers. This includes pointers to the     !
!     computational grid and the momentum matrix.                              !
!   * PISO algorithm implementation: If the simulation is using the PISO       !
!     algorithm for pressure-velocity coupling and is within the PISO loop,    !
!     the subroutine proceeds to update the momentum.                          !
!   * Neighbour summation: Iterates through the faces of the grid to sum       !
!     the contributions from neighboring cells to the momentum.                !
!   * Explicit velocity update: Solves the velocity component (ui) explicitly  !
!     without under-relaxation. This step directly calculates the new velocity !
!     values based on the summed contributions and the momentum equations.     !
!   * Gradient computation: Calculates the gradients of the updated velocity   !
!     field, which is essential for accurately capturing the flow dynamics.    !
!   * Cleanup and performance monitoring: Disconnects any connected real cells !
!     and stops the performance profiler, marking the end of the subroutine.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Var_Type)              :: ui       !! velocity component
  type(Native_Type),   target :: Nat      !! native solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Matrix_Type), pointer :: M
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  real, contiguous,  pointer :: sum_neigh(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Compute_Momentum_Explicit')

  call Work % Connect_Real_Cell(sum_neigh)

  ! Take aliases
  Grid => Flow % pnt_grid
  M    => Nat % M
  b    => Nat % b % val

  ! PISO corrections are executed here
  if (Flow % p_m_coupling == PISO .and. Flow % inside_piso_loop) then

    ! Sum of neighbours
    sum_neigh(:) = 0.0
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      sum_neigh(c1) = sum_neigh(c1) - M % val(M % pos(1,s)) * ui % n(c2)
      sum_neigh(c2) = sum_neigh(c2) - M % val(M % pos(2,s)) * ui % n(c1)
    end do

    ! Solve velocity explicitely (no under relaxation!!)
    do c = Cells_In_Domain()
      ui % n(c) = (sum_neigh(c) + b(c)) / M % val(M % dia(c))
    end do

    call Flow % Grad_Variable(ui)

  end if

  call Work % Disconnect_Real_Cell(sum_neigh)

  call Profiler % Stop('Compute_Momentum_Explicit')

  end subroutine
