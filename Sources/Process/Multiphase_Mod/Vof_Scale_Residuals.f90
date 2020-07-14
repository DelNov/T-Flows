!==============================================================================!
  subroutine Multiphase_Mod_Vof_Scale_Residuals(mult, sol, ini, all_vel)
!------------------------------------------------------------------------------!
!   Scales residuals for PISO algorithm                                        !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod,       only: res => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
  logical                       :: all_vel
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Matrix_Type),pointer :: a
  type(Var_Type),   pointer :: u, v, w
  real, contiguous, pointer :: b(:)
  real                      :: denom, epsloc, simple_tol
  integer                   :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid      => mult % pnt_grid
  flow      => mult % pnt_flow
  a         => sol % a
  b         => sol % b % val
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  epsloc = epsilon(epsloc)

  if (all_vel) then !for velocity
    denom = 0.0
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      denom = denom + abs(a % val(a % dia(c)) *                              &
              norm2((/flow % u % n(c), flow % v % n(c), flow % w % n(c)/)))
    end do

    call Comm_Mod_Global_Sum_Real(denom)

    if (denom <= epsloc) then
      denom = 1.0
    end if
    u % res = u % res_scal / denom
    v % res = v % res_scal / denom
    w % res = w % res_scal / denom

    call Info_Mod_Iter_Fill_At(1, 1, u % name, u % eniter, u % res)
    call Info_Mod_Iter_Fill_At(1, 2, v % name, v % eniter, v % res)
    call Info_Mod_Iter_Fill_At(1, 3, w % name, w % eniter, w % res)
  else ! for temperature
    denom = 0.0
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      denom = denom + abs(a % val(a % dia(c)) * flow % t % n(c))
    end do

    call Comm_Mod_Global_Sum_Real(denom)

    if (denom <= epsloc) then
      denom = 1.0
    end if
    flow % t % res = flow % t % res_scal / denom
  end if

  end subroutine
