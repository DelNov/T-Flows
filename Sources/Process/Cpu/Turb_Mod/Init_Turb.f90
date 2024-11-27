!==============================================================================!
  subroutine Init_Turb(Turb)
!------------------------------------------------------------------------------!
!   Turbulence model initializations (at the beginning of a time step)         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow

  !--------------------------------------------------------!
  !   Start initializations of various turbulence models   !
  !--------------------------------------------------------!

  if(Turb % model .eq. DES_SPALART) then
    call Calculate_Shear_And_Vorticity(Flow)
  end if

  if(Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. LES_DYNAMIC     .or.  &
     Turb % model .eq. LES_WALE) then
    call Calculate_Shear_And_Vorticity(Flow)
    if(Turb % model .eq. LES_DYNAMIC) then
      call Turb % Vis_T_Dynamic()
    end if
    if(Turb % model .eq. LES_WALE) then
      call Turb % Vis_T_Wale()
    end if
    call Turb % Vis_T_Subgrid()
  end if

  if(Turb % model .eq. LES_TVM) then
    call Turb % Vis_T_Tensorial()
  end if

  if(Turb % model .eq. HYBRID_LES_RANS) then
    call Calculate_Shear_And_Vorticity(Flow)
    call Turb % Vis_T_Dynamic()
    call Turb % Vis_T_Hybrid_Les_Rans()
  end if

  if(Turb % model .eq. HYBRID_LES_PRANDTL) then
    call Calculate_Shear_And_Vorticity(Flow)
    call Turb % Vis_T_Hybrid_Les_Prandtl()
  end if

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb % Vis_T_Rsm()
  end if

  end subroutine
