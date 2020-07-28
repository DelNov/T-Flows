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

    flow % m_flux % o(1:) = flow % m_flux % n(1:)  &
                             / flow % density_f(1:)

    ! Update the values at boundaries
    call Update_Boundary_Values(flow, turb, mult)  ! nema veze
    call Multiphase_Mod_Compute_Vof(mult, sol, flow % dt, n)
    call Field_Mod_Body_Forces(flow)
    call Multiphase_Averaging(flow, mult, mult % vof)

    if(mult % track_front) then
    call Surf_Mod_Allocate(mult % surf, flow)
    call Surf_Mod_Place_At_Var_Value(mult % surf,  &
                                     mult % vof,   &
                                     sol,          &
                                     0.5,          &
                                     .true.)  ! don't print messages
    call Surf_Mod_Calculate_Curvatures_From_Elems(mult % surf)
    call Surf_Mod_Compute_Distance_Function_And_Vof(mult % surf,       &
                                                    mult % dist_func,  &
                                                    mult % vof)
    end if
  else
    flow % m_flux % o(1:) = flow % m_flux % n(1:)

  end if

  end subroutine
