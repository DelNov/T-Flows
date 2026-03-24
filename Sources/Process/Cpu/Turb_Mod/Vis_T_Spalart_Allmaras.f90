!==============================================================================!
  subroutine Vis_T_Spalart_Allmaras(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: vis
  integer                   :: c
  real                      :: x_rat, f_v1
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  call Flow % Alias_Momentum(u, v, w)

  !-------------------------!
  !   Turbulent viscosity   !
  !-------------------------+
  do c = Cells_In_Domain_And_Buffers()
    x_rat    = vis % n(c) / (Flow % viscosity(c)/Flow % density(c))
    f_v1     = x_rat**3/(x_rat**3 + Turb % c_v1**3)
    Turb % vis_t(c) = Flow % density(c) * f_v1 * vis % n(c)
    if(Time % Curr_Dt() < 4800)  &
    Turb % vis_t(c) = min(Turb % vis_t(c),2.0*Flow % viscosity(c))
  end do

  !-------------------!
  !   Wall function   !
  !-------------------+
  call Turb % Wall_Function()

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if

  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
