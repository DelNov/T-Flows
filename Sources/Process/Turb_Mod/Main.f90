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

  if(turb % model .eq. K_EPS) then
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

  if(turb % model .eq. K_EPS_ZETA_F .or. &
     turb % model .eq. HYBRID_LES_RANS) then

    call Calculate_Shear_And_Vorticity(flow)

    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % kin, n)
    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % eps, n)

    if(heat_transfer) then
      call Turb_Mod_Calculate_Stress   (turb)
      call Turb_Mod_Calculate_Heat_Flux(turb)
      call Turb_Mod_Compute_Variable(turb, sol, ini, turb % t2, n)
    end if

    call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % zeta, n)

    call Turb_Mod_Vis_T_K_Eps_Zeta_F(turb)

  end if

  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    call Time_And_Length_Scale(grid, turb)

    call Field_Mod_Grad_Variable(flow, flow % u)
    call Field_Mod_Grad_Variable(flow, flow % v)
    call Field_Mod_Grad_Variable(flow, flow % w)

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uu)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vv)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % ww)

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uv)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % uw)
    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % vw)

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Turb_Mod_Compute_F22(turb, sol, ini, turb % f22)
    end if

    call Turb_Mod_Compute_Stress(turb, sol, ini, turb % eps)

    call Turb_Mod_Vis_T_Rsm(turb)

    if(heat_transfer) then
      call Turb_Mod_Calculate_Heat_Flux(turb)
    end if
  end if

  if(turb % model .eq. SPALART_ALLMARAS .or.  &
     turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(flow)
    call Calculate_Vorticity(flow)

    call Turb_Mod_Compute_Variable(turb, sol, ini, turb % vis, n)
    call Turb_Mod_Vis_T_Spalart_Allmaras(turb)
  end if

  end subroutine
