!==============================================================================!
  subroutine Turb_Mod_Main(Turb, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type)     :: Turb
  type(Solver_Type)   :: Sol
  integer, intent(in) :: curr_dt ! current time step
  integer, intent(in) :: ini     ! inner iteration
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(Turb % model .eq. K_EPS) then
    call Calculate_Shear_And_Vorticity(Flow)
    call Time_And_Length_Scale(Grid, Turb)

    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % kin)
    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % eps)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Stress   (Turb)
      call Turb_Mod_Calculate_Heat_Flux(Turb)
      call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % t2)
    end if

    call Turb_Mod_Vis_T_K_Eps(Turb)

  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or. &
     Turb % model .eq. HYBRID_LES_RANS) then

    call Calculate_Shear_And_Vorticity(Flow)

    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % kin)
    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % eps)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Stress   (Turb)
      call Turb_Mod_Calculate_Heat_Flux(Turb)
      call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % t2)
    end if

    call Turb_Mod_Compute_F22(Turb, Sol, curr_dt, ini, Turb % f22)
    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % zeta)

    call Turb_Mod_Vis_T_K_Eps_Zeta_F(Turb)

  end if

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    call Time_And_Length_Scale(Grid, Turb)

    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)

    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % uu)
    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % vv)
    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % ww)

    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % uv)
    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % uw)
    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % vw)

    if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Turb_Mod_Compute_F22(Turb, Sol, curr_dt, ini, Turb % f22)
    end if

    call Turb_Mod_Compute_Stress(Turb, Sol, curr_dt, ini, Turb % eps)

    call Turb_Mod_Vis_T_Rsm(Turb)

    if(Flow % heat_transfer) then
      call Turb_Mod_Calculate_Heat_Flux(Turb)
    end if
  end if

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(Flow)
    call Calculate_Vorticity(Flow)

    call Turb_Mod_Compute_Variable(Turb, Sol, curr_dt, ini, Turb % vis)
    call Turb_Mod_Vis_T_Spalart_Allmaras(Turb)
  end if

  end subroutine
