!==============================================================================!
  subroutine Main_Turb(Turb, Sol)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type)  :: Turb
  type(Solver_Type) :: Sol
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: phi
  integer                   :: sc
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(Turb % model .eq. K_EPS) then

    ! Calculate turbulent scalar fluxes
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      if(Flow % n_scalars > 0) then
        call Turb % Calculate_Stress     ()
        call Turb % Calculate_Scalar_Flux(sc)
      end if
    end do

    call Calculate_Shear_And_Vorticity(Flow)
    call Turb % Time_And_Length_Scale(Grid)

    call Turb % Compute_Variable(Sol, Turb % kin)
    call Turb % Compute_Variable(Sol, Turb % eps)

    if(Flow % heat_transfer) then
      call Turb % Calculate_Stress   ()
      call Turb % Calculate_Heat_Flux()
      call Turb % Compute_Variable(Sol, Turb % t2)
    end if

    call Turb % Vis_T_K_Eps()

  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or. &
     Turb % model .eq. HYBRID_LES_RANS) then

    ! Calculate turbulent scalar fluxes
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      if(Flow % n_scalars > 0) then
        call Turb % Calculate_Stress     ()
        call Turb % Calculate_Scalar_Flux(sc)
      end if
    end do

    call Calculate_Shear_And_Vorticity(Flow)

    call Turb % Compute_Variable(Sol, Turb % kin)
    call Turb % Compute_Variable(Sol, Turb % eps)

    if(Flow % heat_transfer) then
      call Turb % Calculate_Stress   ()
      call Turb % Calculate_Heat_Flux()
      call Turb % Compute_Variable(Sol, Turb % t2)
    end if

    call Turb % Compute_F22(Sol, Turb % f22)
    call Turb % Compute_Variable(Sol, Turb % zeta)

    ! For some cases, it is beneficial to start simulations with
    ! turbulent viscosity computed with k-eps.  Particularly for
    ! cases with mild pressure drops such as channel, pipe flows
    if(Time % Curr_Dt() < 10) then
      call Turb % Vis_T_K_Eps()
    else
      call Turb % Vis_T_K_Eps_Zeta_F()
    end if

  end if

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    call Turb % Time_And_Length_Scale(Grid)

    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)

    call Turb % Compute_Stress(Sol, Turb % uu)
    call Turb % Compute_Stress(Sol, Turb % vv)
    call Turb % Compute_Stress(Sol, Turb % ww)

    call Turb % Compute_Stress(Sol, Turb % uv)
    call Turb % Compute_Stress(Sol, Turb % uw)
    call Turb % Compute_Stress(Sol, Turb % vw)

    if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Turb % Compute_F22(Sol, Turb % f22)
    end if

    call Turb % Compute_Stress(Sol, Turb % eps)

    call Turb % Vis_T_Rsm()

    if(Flow % heat_transfer) then
      call Turb % Calculate_Heat_Flux()
    end if
  end if

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(Flow)

    call Turb % Compute_Variable(Sol, Turb % vis)
    call Turb % Vis_T_Spalart_Allmaras()
  end if

  end subroutine
