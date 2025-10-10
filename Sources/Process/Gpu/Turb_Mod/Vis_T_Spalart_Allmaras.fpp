!==============================================================================!
  subroutine Vis_T_Spalart_Allmaras(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Grid_Type)          :: Grid
  type(Field_Type)         :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),   pointer :: u, v, w, t
  type(Var_Type),   pointer :: vis
  integer                   :: c
  real                      :: x_rat, f_v1
!==============================================================================!

  ! Take aliases
  vis  => Turb % vis
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  !-----------------------------------!
  !   Calculate turbulent viscosity   !
  !-----------------------------------!
  if(Turb % model .eq. DES_SPALART) then
    do c = Cells_In_Domain_And_Buffers()
      x_rat    = vis % n(c) / (Flow % viscosity(c)/Flow % density(c))
      f_v1     = x_rat**3/(x_rat**3 + Turb % c_v1**3)
      Turb % vis_t(c) = Flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  if(Turb % model .eq. SPALART_ALLMARAS) then
    do c = Cells_In_Domain_And_Buffers()
      x_rat = vis % n(c) / (Flow % viscosity(c)/Flow % density(c))
      f_v1  = x_rat**3/(x_rat**3 + Turb % c_v1**3)
      Turb % vis_t(c) = Flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  !-------------------!
  !   Wall function   !
  !-------------------!
  call Turb % Wall_Function(Grid, Flow)

  ! Refresh buffers
  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
