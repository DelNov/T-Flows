!==============================================================================!
  subroutine Surface_Tension_Force(Vof, i)
!------------------------------------------------------------------------------!
!   Computes Surface tension, forces for momentum                              !
!                                                                              !
!   Recall how momentum equations looks in Compute_Momentum procedure:         !
!                                                                              !
!     /             /              /               /             /             !
!    |     du      |              |               |             |              !
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | GRAD p dV + | f dV         !
!    |     dt      |              |               |             |              !
!   /             /              /               /             /               !
!                                                                              !
!   Here, f will be the surface tension force.  If we use material             !
!   derivative for the sake of shorteness, the equation will read:             !
!                                                                              !
!     /             /               /             /                            !
!    |     Du      |               |             |                             !
!    | rho -- dV = | mu DIV u dS - | GRAD p dV + | f dV                        !
!    |     Dt      |               |             |                             !
!   /             /               /             /                              !
!                                                                              !
!   and with definition of surface tension force:                              !
!                                                                              !
!     /             /               /                        /                 !
!    |     Du      |               |                         |                 !
!    | rho -- dV = | mu DIV u dS - | GRAD p dV + sigma kappa | GRAD c dV       !
!    |     Dt      |               |                         |                 !
!   /             /               /                         /                  !
!                                                                              !
!   In order to balance pressure with surface tension force, last two terms    !
!   should always be treated in exactly the same manner!  One could even       !
!   write the equation like this:                                              !
!                                                                              !
!     /             /               /                                          !
!    |     Du      |               |                                           !
!    | rho -- dV = | mu DIV u dS - | GRAD ( p + sigma kappa c ) dV             !
!    |     Dt      |               |                                           !
!   /             /               /       |                   |                !
!                                         |<---- together --->|                !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer, intent(in)     :: i
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: col
  real, contiguous,  pointer :: surf_fx(:), surf_fy(:), surf_fz(:)
  integer                    :: c
!==============================================================================!

  ! Get out of here if surface tension is neglected
  if(Vof % surface_tension < TINY) return

  ! Take aliases
  Flow    => Vof % pnt_flow
  Grid    => Vof % pnt_grid
  col     => Vof % fun
  surf_fx => Vof % surf_fx
  surf_fy => Vof % surf_fy
  surf_fz => Vof % surf_fz
  ! Don't do what you once did in the line which follows this comment,
  ! it is plain silly as the smoothed (convoluted) variant of the vof
  ! function is used just for estimation of normals and curvatures.
  ! col    => Vof % smooth

  ! Units here:
  ! N/m * 1/m * 1/m = N / m^3
  select case(i)
    case(1)
      do c = 1, Grid % n_cells
        surf_fx(c) = Vof % surface_tension  &
                   * Vof % curv(c)          &
                   * col % x(c)
       end do
    case(2)
      do c = 1, Grid % n_cells
        surf_fy(c) = Vof % surface_tension  &
                   * Vof % curv(c)          &
                   * col % y(c)
      end do
    case(3)
      do c = 1, Grid % n_cells
        surf_fz(c) = Vof % surface_tension  &
                   * Vof % curv(c)          &
                   * col % z(c)
      end do

  end select

  end subroutine
