!==============================================================================!
  subroutine Multiphase_Mod_Main(mult, flow, turb, sol, n)
!------------------------------------------------------------------------------!
!   Initialize inteface tracking simulations                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Solver_Type)             :: sol
  integer, intent(in)           :: n     ! time step
!==============================================================================!

  ! Volume of Fluid
  if(mult % model .eq. VOLUME_OF_FLUID) then

    flow % v_flux % n(1:) = flow % m_flux % n(1:) / flow % density_f(1:)

    ! Front tracking is not engaged
    if(.not. mult % track_front) then
      call Update_Boundary_Values(flow, turb, mult)
      call Multiphase_Mod_Compute_Vof(mult, sol, flow % dt, n)
      call Field_Mod_Body_Forces(flow)
      call Multiphase_Averaging(flow, mult, mult % vof)

    ! Front tracking is engaged
    else
      call Update_Boundary_Values(flow, turb, mult)
    ! call Surf_Mod_Advance_Front(mult, flow % dt, n) not implemented yet
      call Field_Mod_Body_Forces(flow)
    end if

  else
    flow % v_flux % n(1:) = flow % m_flux % n(1:)

  end if

  end subroutine
