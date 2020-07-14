!============================================================================!
  subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, sol, ui, i)
!----------------------------------------------------------------------------!
!   Computes Surface tension, Gravity and phase change sources for Momentum  !
!   Equation if a two-phase flow calculation is performed. Additionally and  !
!   for the moment, PISO calculations are run here                           !
!----------------------------------------------------------------------------!
!--------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: neigh => r_cell_12,   &
                      res   => r_cell_13
!----------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]--------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  type(Var_Type),        target :: ui
  integer                       :: i
!-----------------------------------[Locals]---------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  type(Face_Type),  pointer :: m_flux
  type(Matrix_Type),pointer :: a
  real, contiguous, pointer :: b(:)
  integer                   :: s, c, c1, c2, nt, ni
  real                      :: dotprod, epsloc, fs
  real,             pointer :: si(:), sj(:), sk(:)
  real                      :: corr_x, corr_y, corr_z
  real                      :: u_f, v_f, w_f
!============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => mult % pnt_grid
  vof    => mult % vof
  m_flux => flow % m_flux
  si     => grid % sx
  sj     => grid % sy
  sk     => grid % sz
  a      => sol % a
  b      => sol % b % val

  epsloc = epsilon(epsloc)

  ! Surface tension contribution
  if (mult % surface_tension > TINY) then

    select case(i)
      case(1)
        do c = 1, grid % n_cells
          b(c) = b(c) + mult % surface_tension  &
                      * mult % curv(c)          &
                      * vof % x(c)              &
                      * grid % vol(c)
         end do
      case(2)
        do c = 1, grid % n_cells
          b(c) = b(c) + mult % surface_tension  &
                      * mult % curv(c)          &
                      * vof % y(c)              &
                      * grid % vol(c)
        end do
      case(3)
        do c = 1, grid % n_cells
          b(c) = b(c) + mult % surface_tension  &
                      * mult % curv(c)          &
                      * vof % z(c)              &
                      * grid % vol(c)
        end do

    end select

  end if

  ! Body force contribution (gravity)
  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    select case(i)

      case(1)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fx(c)
        end do

      case(2)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fy(c)
        end do

      case(3)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fz(c)
        end do

    end select

  end if

  ! Momentum variables for pressure correction
  ! This is here because they need to be collected before
  ! u, v, w are calculated

  if (flow % temp_corr) then
    ! Guessed face velocity
    if (i == 1) then
      do s = grid % n_bnd_faces + 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        fs = grid % f(s)
        u_f = fs * flow % u % n(c1) + (1.0 - fs) * flow % u % n(c2)
        v_f = fs * flow % v % n(c1) + (1.0 - fs) * flow % v % n(c2)
        w_f = fs * flow % w % n(c1) + (1.0 - fs) * flow % w % n(c2)
        m_flux % avg(s) = ( u_f * grid % sx(s)     &
                          + v_f * grid % sy(s)     &
                          + w_f * grid % sz(s) )
      end do

      do s = grid % n_bnd_faces + 1, grid % n_faces
        m_flux % star(s) = m_flux % n(s) / flow % density_f(s)
      end do
    end if
  end if

  end subroutine
