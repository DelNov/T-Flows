!============================================================================!
  subroutine Surface_Tension_Force(Vof, fi, i)
!----------------------------------------------------------------------------!
!   Computes Surface tension, forces for momentum                            !
!----------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]--------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: fi(:)
  integer, intent(in)     :: i
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
  flow    => Vof % pnt_flow
  grid    => Vof % pnt_grid
  ! col     => Vof % smooth
  col     => Vof % fun
  surf_fx => Vof % surf_fx
  surf_fy => Vof % surf_fy
  surf_fz => Vof % surf_fz
  v_flux  => flow % v_flux

  ! Surface tension contribution
  if(Vof % surface_tension > TINY) then

    ! Units here:
    ! N/m * 1/m * 1/m * m^3 = N
    select case(i)
      case(1)
        do c = 1, grid % n_cells
          surf_fx(c) = Vof % surface_tension  &
                     * Vof % curv(c)          &
                     * col % x(c)             &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fx(c)
         end do
      case(2)
        do c = 1, grid % n_cells
          surf_fy(c) = Vof % surface_tension  &
                     * Vof % curv(c)          &
                     * col % y(c)             &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fy(c)
        end do
      case(3)
        do c = 1, grid % n_cells
          surf_fz(c) = Vof % surface_tension  &
                     * Vof % curv(c)          &
                     * col % z(c)             &
                     * grid % vol(c)
          fi(c) = fi(c) + surf_fz(c)
        end do

    end select

  end if

  end subroutine
