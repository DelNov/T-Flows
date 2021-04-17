!============================================================================!
  subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, sol, i)
!----------------------------------------------------------------------------!
!   Computes Surface tension, Gravity and phase change sources for Momentum  !
!   Equation if a two-phase flow calculation is performed. Additionally and  !
!   for the moment, PISO calculations are run here                           !
!----------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]--------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: i
!-----------------------------------[Locals]---------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: col
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: m
  real, contiguous,  pointer :: b(:)
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
  m       => sol % m
  b       => sol % b % val

  ! Surface tension contribution
  if(mult % surface_tension > TINY) then

    select case(i)
      case(1)
        do c = 1, grid % n_cells
          surf_fx(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % x(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fx(c)
         end do
      case(2)
        do c = 1, grid % n_cells
          surf_fy(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % y(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fy(c)
        end do
      case(3)
        do c = 1, grid % n_cells
          surf_fz(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * col % z(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fz(c)
        end do

    end select

  end if

  ! Momentum variables for pressure correction
  ! This is here because they need to be collected
  ! before u, v, w are calculated
  ! Bojan: This is very strange, works for one velocity component only!
  if(flow % temp_corr) then

    ! Guessed face velocity
    if(i == 1) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        if(c2 > 0) then
          fs = grid % f(s)
          u_f = fs * flow % u % n(c1) + (1.0 - fs) * flow % u % n(c2)
          v_f = fs * flow % v % n(c1) + (1.0 - fs) * flow % v % n(c2)
          w_f = fs * flow % w % n(c1) + (1.0 - fs) * flow % w % n(c2)
          v_flux % avg(s) = (  u_f * grid % sx(s)     &
                             + v_f * grid % sy(s)     &
                             + w_f * grid % sz(s) )
        end if
      end do

      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        if(c2 > 0) v_flux % star(s) = v_flux % n(s)
      end do
    end if
  end if

  end subroutine
