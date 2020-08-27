!==============================================================================!
  subroutine Piso_Algorithm(flow, turb, mult, sol, mass_res, ini)
!------------------------------------------------------------------------------!
!   PISO algorithm                                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Turb_Mod
  use Multiphase_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Calling]------------------------------------!
  real :: Correct_Velocity
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini       ! current inner iteration
  real                          :: mass_res
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: a
  type(Grid_Type),   pointer :: grid
  real, contiguous,  pointer :: b(:)
  integer                    :: corr_steps
!==============================================================================!

  grid => flow % pnt_grid
  call Solver_Mod_Alias_System(sol, a, b)

  if (flow % p_m_coupling == PISO) then

    flow % piso_status = .true.
    do corr_steps = 1, flow % n_piso_corrections
      flow % i_corr = corr_steps
      call Field_Mod_Grad_Pressure(flow, flow % p)

      ! All velocity components one after another
      call Compute_Momentum(flow, turb, mult, sol, flow % dt, ini)

      ! Refresh buffers for a % sav before discretizing for pressure
      ! (Can this call be somewhere in Compute Pressure?)
      call Grid_Mod_Exchange_Cells_Real(grid, sol % a % sav)

      call Balance_Volume(flow, mult)
      call Compute_Pressure(flow, mult, sol, flow % dt, ini)
      call Multiphase_Averaging(flow, mult, flow % p)

      mass_res = Correct_Velocity(flow, mult, sol, flow % dt, ini)
      call Multiphase_Averaging(flow, mult, flow % u)
      call Multiphase_Averaging(flow, mult, flow % v)
      call Multiphase_Averaging(flow, mult, flow % w)
    end do
    flow % piso_status = .false.
    call Multiphase_Mod_Vof_Scale_Residuals(mult, sol, ini, .true.)
    flow % i_corr = 1
  end if

  end subroutine
