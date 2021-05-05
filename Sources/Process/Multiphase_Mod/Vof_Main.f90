!==============================================================================!
  subroutine Multiphase_Mod_Vof_Main(mult, flow, turb, Sol, n)
!------------------------------------------------------------------------------!
!   Initialize inteface tracking simulations                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Solver_Type)             :: Sol
  integer, intent(in)           :: n     ! time step
!==============================================================================!

  flow % v_flux % o(1:) = flow % v_flux % n(1:)

  ! Volume of Fluid
  if(mult % model .eq. VOLUME_OF_FLUID) then
    if(mult % track_front) then
      call Front_Mod_Place_At_Var_Value(mult % front,  &
                                        mult % vof,    &
                                        Sol,           &
                                        0.5,           &
                                        .true.)  ! don't print messages
      call Front_Mod_Calculate_Curvatures_From_Elems(mult % front)
      call Front_Mod_Print_Statistics               (mult % front)
!f_vs_s      call Surf_Mod_Place_At_Var_Value(mult % surf,  &
!f_vs_s                                       mult % vof,   &
!f_vs_s                                       Sol,          &
!f_vs_s                                       0.5,          &
!f_vs_s                                       .true.)  ! don't print messages
!f_vs_s      call Surf_Mod_Calculate_Curvatures_From_Elems(mult % surf)
!f_vs_s      call Surf_Mod_Print_Statistics               (mult % surf)
    end if

    ! Commands to advance vof
    call Update_Boundary_Values(flow, turb, mult, 'MULTIPHASE')
    call Multiphase_Mod_Vof_Compute(mult, Sol, flow % dt, n)

  end if

  end subroutine
