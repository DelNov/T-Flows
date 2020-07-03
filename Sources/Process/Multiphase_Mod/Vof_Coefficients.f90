!==============================================================================!
  subroutine Multiphase_Mod_Vof_Coefficients(mult, a, b, dt, beta_f)
!------------------------------------------------------------------------------!
!   Computes matrix coefficients for volume fraction equation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Matrix_Type),     target :: a
  real,                  target :: b(:)
  real                          :: dt
  real                          :: beta_f(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: m_flux
  type(Var_Type),   pointer :: vof
  integer                   :: c, c1, c2, s
  real                      :: upwd1, upwd2, upwd3, a0
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  m_flux => flow % m_flux
  vof    => mult % vof

  ! Initialize matrix and right hand side
  b       = 0.0
  a % val = 0.0

  !-------------------------!
  !   Matrix Coefficients   !
  !-------------------------!
  if (vof % adv_scheme .eq. UPWIND) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      upwd1 = 0.5 * max( m_flux % n(s) / flow % density_f(s), 0.0)
      upwd2 = 0.5 * max(-m_flux % n(s) / flow % density_f(s), 0.0)

      a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1
      b(c1) = b(c1) - ( upwd1 * vof % o(c1) -  upwd2 * vof % o(c2) ) 

      if (c2 > 0) then
        a % val(a % pos(1,s)) =  - upwd2

        a % val(a % dia(c2)) = a % val(a % dia(c2)) + upwd2
        b(c2) = b(c2) - ( upwd2 * vof % o(c2) - upwd1 * vof % o(c1))
        a % val(a % pos(2,s)) =  - upwd1
      else
        if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          b(c1) = b(c1) - m_flux % n(s) / flow % density_f(s) * vof % n(c2)
        else if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          a % val(a % dia(c1)) = a % val(a % dia(c1))          &
                               + m_flux % n(s) / flow % density_f(s)
        end if
      end if

    end do

  else if (vof % adv_scheme .eq. CICSAM .or. &
           vof % adv_scheme .eq. STACS) then

    do s = 1, grid % n_faces

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      upwd1 = (0.5 - beta_f(s)) * max(-m_flux % n(s)                  &
                                     / flow % density_f(s), 0.0)      &
             - 0.5 * beta_f(s) * m_flux % n(s) / flow % density_f(s)
      upwd2 = (0.5 - beta_f(s)) * max(m_flux % n(s)                   &
                                     / flow % density_f(s), 0.0)      &
             + 0.5 * beta_f(s) * m_flux % n(s) / flow % density_f(s)
      upwd3 = 0.5 * m_flux % n(s) / flow % density_f(s)

      if (c2 > 0) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1 + upwd3
        a % val(a % pos(1,s)) = -upwd1
        b(c1) = b(c1) - (upwd1 + upwd3) * vof % o(c1) + upwd1 * vof % o(c2)

        a % val(a % dia(c2)) = a % val(a % dia(c2)) + upwd2 - upwd3
        a % val(a % pos(2,s)) = -upwd2 
        b(c2) = b(c2) - (upwd2 - upwd3) * vof % o(c2) + upwd2 * vof % o(c1)
      else
        b(c1) = b(c1) - m_flux % n(s) / flow % density_f(s) * vof % n(c2)
      end if
    end do

  end if

  !------------------------------------------------------------!
  !                Calculate Main coefficient and              !
  !   Calculate Source term (part 2 : Temporal contribution)   !
  !------------------------------------------------------------!

  ! Two time levels; Linear interpolation
  if(vof % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * vof % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(vof % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0 * a0 * vof % o(c) - 0.5 * a0 * vof % oo(c)
    end do
  end if

  end subroutine
