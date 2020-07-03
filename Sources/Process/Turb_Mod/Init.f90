!==============================================================================!
  subroutine Turb_Mod_Init(turb)
!------------------------------------------------------------------------------!
!   Turbulence model initializations (at the beginning of a time step)         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: flow
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow

  !--------------------------------------------------------!
  !   Start initializations of various turbulence models   !
  !--------------------------------------------------------!

  if(turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(flow)
    call Calculate_Vorticity (flow)
  end if

  if(turb % model .eq. LES_SMAGORINSKY .or.  &
     turb % model .eq. LES_DYNAMIC     .or.  &
     turb % model .eq. LES_WALE) then
    call Calculate_Shear_And_Vorticity(flow)
    if(turb % model .eq. LES_DYNAMIC) then
      call Turb_Mod_Vis_T_Dynamic(turb)
    end if
    if(turb % model .eq. LES_WALE) then
      call Turb_Mod_Vis_T_Wale(turb)
    end if
    call Turb_Mod_Vis_T_Smagorinsky(turb)
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    call Calculate_Shear_And_Vorticity(flow)
    call Turb_Mod_Vis_T_Dynamic(turb)
    call Turb_Mod_Vis_T_Hybrid_Les_Rans(turb)
  end if

  if(turb % model .eq. HYBRID_LES_PRANDTL) then
    call Calculate_Shear_And_Vorticity(flow)
    call Turb_Mod_Vis_T_Hybrid_Les_Prandtl(turb)
  end if

  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb_Mod_Vis_T_Rsm(turb)
  end if

  end subroutine
