!==============================================================================!
  subroutine Turb_Mod_Main(turb, sol, n, ini)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type)   :: turb
  type(Solver_Type) :: sol
  integer           :: n     ! time step
  integer           :: ini   ! inner iteration
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(turbulence_model .eq. K_EPS) then

    ! Update the values at boundaries
    call Update_Boundary_Values(flow, turb)

    call Calculate_Shear_And_Vorticity(flow)

    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % kin, n)
    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % eps, n)

    if(heat_transfer) then
      call Turb_Mod_Calculate_Stress   (turb)
      call Turb_Mod_Calculate_Heat_Flux(turb)
      call Turb_Mod_Compute_Variable(turb, sol, ini, turb % t2, n)
    end if

    call Turb_Mod_Vis_T_K_Eps(turb)

  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .or. &
     turbulence_model .eq. HYBRID_LES_RANS) then
    call Calculate_Shear_And_Vorticity(flow)

    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % kin, n)
    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % eps, n)

    if(heat_transfer) then
      call Turb_Mod_Calculate_Stress   (turb)
      call Turb_Mod_Calculate_Heat_Flux(turb)
      call Turb_Mod_Compute_Variable(turb, sol, ini, turb % t2, n)
    end if

    call Update_Boundary_Values(flow, turb)

    call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % zeta, n)

    call Turb_Mod_Vis_T_K_Eps_Zeta_F(turb)

  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Update the values at boundaries
    call Update_Boundary_Values(flow, turb)

    call Time_And_Length_Scale(grid, turb)

    call Field_Mod_Grad_Variable(flow, flow % u)
    call Field_Mod_Grad_Variable(flow, flow % v)
    call Field_Mod_Grad_Variable(flow, flow % w)

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uu, n)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vv, n)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % ww, n)

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uv, n)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uw, n)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vw, n)

    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
    end if

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % eps, n)

    call Turb_Mod_Vis_T_Rsm(turb)

    if(heat_transfer) then
      call Turb_Mod_Calculate_Heat_Flux(turb)
    end if
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(flow)
    call Calculate_Vorticity(flow)

    ! Update the values at boundaries
    call Update_Boundary_Values(flow, turb)

    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % vis, n)
    call Turb_Mod_Vis_T_Spalart_Allmaras(turb)
  end if

  end subroutine
