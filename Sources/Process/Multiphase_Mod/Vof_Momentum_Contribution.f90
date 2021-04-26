!============================================================================!
  subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, fi, i)
!----------------------------------------------------------------------------!
!   Computes Surface tension, Gravity and phase change sources for Momentum  !
!   Equation if a two-phase flow calculation is performed. Additionally and  !
!   for the moment, PISO calculations are run here                           !
!----------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]--------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: fi(:)
  integer, intent(in)           :: i
!-----------------------------------[Locals]---------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: col
  type(Face_Type),   pointer :: v_flux
  real, contiguous,  pointer :: surf_fx(:), surf_fy(:), surf_fz(:)
  integer                    :: s, c, c1, c2
  real                       :: fs
  real                       :: u_f, v_f, w_f
!============================================================================!

  ! Take aliases
  flow    => mult % pnt_flow
  grid    => mult % pnt_grid
  ! col     => mult % smooth
  col     => mult % vof
  surf_fx => mult % surf_fx
  surf_fy => mult % surf_fy
  surf_fz => mult % surf_fz
  v_flux  => flow % v_flux

  ! Surface tension contribution
  if(mult % surface_tension > TINY) then

    select case(i)
      case(1)
        do c = 1, grid % n_cells
          surf_fx(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % x(c)              &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fx(c)
         end do
      case(2)
        do c = 1, grid % n_cells
          surf_fy(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % y(c)              &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fy(c)
        end do
      case(3)
        do c = 1, grid % n_cells
          surf_fz(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % z(c)              &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fz(c)
        end do

    end select

  end if

  end subroutine
