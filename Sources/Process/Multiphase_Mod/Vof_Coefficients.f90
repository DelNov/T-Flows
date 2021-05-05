!==============================================================================!
  subroutine Multiphase_Mod_Vof_Coefficients(mult, a, b, dt)
!------------------------------------------------------------------------------!
!   Computes matrix coefficients for volume fraction equation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Matrix_Type),     target :: a
  real,                  target :: b(:)
  real                          :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: v_flux
  type(Var_Type),   pointer :: vof
  type(Front_Type), pointer :: front
  real, contiguous, pointer :: beta_f(:)
  integer                   :: c, c1, c2, s, l, g, e
  real                      :: upwd1, upwd2, upwd3, a0
  real                      :: vof_int, x_int, y_int, z_int, w_c1
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  v_flux => flow % v_flux
  vof    => mult % vof
  beta_f => mult % beta_f
  front  => mult % front

  ! Distinguish between liquid and vapor
  l = 1; g = 2
  if(mult % phase_dens(g) > mult % phase_dens(l)) then
    l = 2; g = 1
  end if

  ! Mark cells at surface
  if(flow % mass_transfer) then
    mult % cell_at_elem(:) = 0  ! not at surface
    do e = 1, front % n_elems
      mult % cell_at_elem(front % elem(e) % cell) = e
    end do
  end if

  ! Initialize matrix and right hand side
  b       = 0.0
  A % val = 0.0

  !-------------------------!
  !   Matrix Coefficients   !
  !-------------------------!
  if(vof % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          b(c1) = b(c1) - v_flux % n(s) * vof % n(c2)
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + v_flux % n(s)
        end if
      end if

    end do

    ! Interior faces
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 > 0) then

        upwd1 = 0.5 * max( v_flux % n(s), 0.0)
        upwd2 = 0.5 * max(-v_flux % n(s), 0.0)

        A % val(A % dia(c1)) = A % val(A % dia(c1)) + upwd1
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + upwd2

        A % val(A % pos(1,s)) =  - upwd2
        A % val(A % pos(2,s)) =  - upwd1

        b(c1) = b(c1) - ( upwd1 * vof % o(c1) - upwd2 * vof % o(c2) )
        b(c2) = b(c2) - ( upwd2 * vof % o(c2) - upwd1 * vof % o(c1) )

      end if
    end do

  else if(vof % adv_scheme .eq. CICSAM .or. &
          vof % adv_scheme .eq. STACS) then

    ! At boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then

        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + v_flux % n(s)
        else
          b(c1) = b(c1) - v_flux % n(s) * vof % n(c2)
        end if

      end if
    end do

    ! Interior faces
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 > 0) then

        upwd1 = (0.5 - beta_f(s)) * max(-v_flux % n(s), 0.0)  &
               - 0.5 * beta_f(s)  * v_flux % n(s)
        upwd2 = (0.5 - beta_f(s)) * max(+v_flux % n(s), 0.0)  &
               + 0.5 * beta_f(s)  * v_flux % n(s)
        upwd3 = 0.5 * v_flux % n(s)

        A % val(A % dia(c1)) = A % val(A % dia(c1)) + upwd1 + upwd3
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + upwd2 - upwd3

        A % val(A % pos(1,s)) = - upwd1
        A % val(A % pos(2,s)) = - upwd2

        b(c1) = b(c1) - (upwd1 + upwd3) * vof % o(c1) + upwd1 * vof % o(c2)
        b(c2) = b(c2) - (upwd2 - upwd3) * vof % o(c2) + upwd2 * vof % o(c1)

      end if
    end do

  end if

  !------------------------------------------------!
  !   Calculate main coefficient and source term   !
  !        (Part 2 : temporal contribution)        !
  !------------------------------------------------!

  ! Two time levels; linear interpolation
  if(vof % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + a0
      b(c) = b(c) + a0 * vof % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(vof % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0 * a0 * vof % o(c) - 0.5 * a0 * vof % oo(c)
    end do
  end if

  ! Phase change
  if(flow % mass_transfer) then
    call Multiphase_Mod_Vof_Mass_Transfer_Vof_Source(mult, b)
  end if

  end subroutine
