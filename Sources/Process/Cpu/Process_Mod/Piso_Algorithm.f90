!==============================================================================!
  subroutine Piso_Algorithm(Process, Flow, Turb, Vof, Por, Sol)
!------------------------------------------------------------------------------!
!>  The subroutine implements the PISO algorithm. This algorithm is designed
!>  for solving fluid flow problems, particularly focusing on pressure-velocity
!>  coupling in CFD.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Sets up the necessary variables and pointers, such as    !
!     grid and flow field, for the PISO algorithm execution. It also           !
!     sets aliases for momentum variables (u, v, w) for shorter referencing.   !
!   * Iterative process: Executes the core steps of the PISO algorithm within  !
!     a loop, iterating through a predefined number of correction steps.       !
!     - Momentum computation: Calculates the momentum conservation equations   !
!       for the fluid flow.                                                    !
!     - Pressure computation: Solves the pressure equation using methods       !
!       specific to the SIMPLE algorithm, a precursor to the PISO algorithm.   !
!     - Velocity correction: Adjusts the velocity fields based on the computed !
!       pressure fields to ensure non-divergence of flow and consistent        !
!       pressure-velocity coupling.                                            !
!   * Post-iteration processing: After the iterative process, updates          !
!     iteration information and resets the correction step counter.            !
!   * User-defined functions: Incorporates user-defined functions to allow for !
!     specific customization and behaviors within the algorithm execution.     !
!   * Performance monitoring: Monitors the algorithm's performance throughout  !
!     its execution, contributing to analysis and optimization of the code.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Turb_Type),     target :: Turb     !! turbulence object
  type(Vof_Type),      target :: Vof      !! VOF object
  type(Porosity_Type), target :: Por      !! porosity object
  type(Solver_Type),   target :: Sol      !! solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: corr_steps
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  if (Flow % p_m_coupling == PISO) then

    Flow % inside_piso_loop = .true.
    do corr_steps = 1, Flow % n_piso_corrections
      Flow % i_corr = corr_steps

      call Process % Compute_Momentum(Flow, Turb, Vof, Por, Sol)
      call Process % Compute_Pressure(Flow,       Vof,      Sol)
      call Process % Correct_Velocity(Flow,       Vof,      Sol)
    end do

    Flow % inside_piso_loop = .false.
    call Info % Iter_Fill_At(1, 1, u % name, u % res, u % niter)
    call Info % Iter_Fill_At(1, 2, v % name, v % res, v % niter)
    call Info % Iter_Fill_At(1, 3, w % name, w % res, w % niter)
    Flow % i_corr = 1

  end if

  end subroutine
