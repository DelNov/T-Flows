!==============================================================================!
  subroutine Turb_Mod_Main(turb, Nat, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type)     :: turb
  type(Native_Type)   :: Nat
  integer, intent(in) :: curr_dt ! current time step
  integer, intent(in) :: ini     ! inner iteration
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(turb % model .eq. K_EPS) then
    call Calculate_Shear_And_Vorticity(Flow)

    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % kin)
    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % eps)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Stress   (turb)
      call Turb_Mod_Calculate_Heat_Flux(turb)
      call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % t2)
    end if

    call Turb_Mod_Vis_T_K_Eps(turb)

  end if

  if(turb % model .eq. K_EPS_ZETA_F .or. &
     turb % model .eq. HYBRID_LES_RANS) then

    call Calculate_Shear_And_Vorticity(Flow)

    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % kin)
    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % eps)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Stress   (turb)
      call Turb_Mod_Calculate_Heat_Flux(turb)
      call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % t2)
    end if

    call Turb_Mod_Compute_F22(turb, Nat, curr_dt, ini, turb % f22)
    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % zeta)

    call Turb_Mod_Vis_T_K_Eps_Zeta_F(turb)

  end if

  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    call Time_And_Length_Scale(Grid, turb)

    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)

    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % uu)
    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % vv)
    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % ww)

    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % uv)
    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % uw)
    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % vw)

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Turb_Mod_Compute_F22(turb, Nat, curr_dt, ini, turb % f22)
    end if

    call Turb_Mod_Compute_Stress(turb, Nat, curr_dt, ini, turb % eps)

    call Turb_Mod_Vis_T_Rsm(turb)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Heat_Flux(turb)
    end if
  end if

  if(turb % model .eq. SPALART_ALLMARAS .or.  &
     turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(Flow)
    call Calculate_Vorticity(Flow)

    call Turb_Mod_Compute_Variable(turb, Nat, curr_dt, ini, turb % vis)
    call Turb_Mod_Vis_T_Spalart_Allmaras(turb)
  end if

  end subroutine
