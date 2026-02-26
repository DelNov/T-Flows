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

  ! Calculate turbulent scalar fluxes
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    if(Flow % n_scalars > 0) then
      call Turb % Calculate_Stress     ()
      call Turb % Calculate_Scalar_Flux(sc)
    end if
  end do

  call Calculate_Shear_And_Vorticity(Flow)

  if(Turb % model .eq. K_EPS) then

    call Turb % Time_And_Length_Scale(Grid)

    call Turb % Compute_Variable(Sol, Turb % kin)
    call Turb % Compute_Variable(Sol, Turb % eps)

    call Turb % Vis_T_K_Eps()

  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or. &
     Turb % model .eq. HYBRID_LES_RANS) then

    call Turb % Compute_Variable(Sol, Turb % kin)
    call Turb % Compute_Variable(Sol, Turb % eps)

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

  if(Turb % model .eq. K_OMEGA_SST) then

    call Turb % Time_And_Length_Scale(Grid)

    call Turb % Update_K_Omega_Sst_Fields()  

    call Turb % Compute_Variable(Sol, Turb % kin)
    call Turb % Compute_Variable(Sol, Turb % omega)

    call Turb % Update_K_Omega_Sst_Fields()
    call Turb % Vis_T_K_Omega_Sst()

  end if


  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then

    call Turb % Compute_Variable(Sol, Turb % vis)
    call Turb % Vis_T_Spalart_Allmaras()

  end if

  if(Flow % heat_transfer.and.Turb % model.ne.NO_TURBULENCE_MODEL) then

    call Turb % Calculate_Stress   ()
    call Turb % Calculate_Heat_Flux()
    if(Flow % buoyancy .eq. THERMALLY_DRIVEN.or.  &
       Turb % heat_flux_model .eq. AFM)           &     
    call Turb % Compute_Variable(Sol, Turb % t2)

  end if

  end subroutine
