!==============================================================================!
  subroutine Discretize(Vof, A, b, dt)
!------------------------------------------------------------------------------!
!   Computes matrix coefficients for volume fraction equation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Matrix_Type), target :: A
  real,              target :: b(:)
  real                      :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Face_Type),  pointer :: v_flux
  type(Var_Type),   pointer :: fun
  type(Front_Type), pointer :: Front
  real, contiguous, pointer :: beta_f(:)
  integer                   :: c, c1, c2, s
  real                      :: upwd1, upwd2, upwd3, a0
!==============================================================================!

  ! Take aliases
  Flow   => Vof % pnt_flow
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  fun    => Vof % fun
  beta_f => Vof % beta_f
  Front  => Vof % Front

  ! Initialize matrix and right hand side
  b       = 0.0
  A % val = 0.0

  !-------------------------!
  !   Matrix Coefficients   !
  !-------------------------!
  if(fun % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW) then
          b(c1) = b(c1) - v_flux % n(s) * fun % n(c2)
        else if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + v_flux % n(s)
        end if
      end if

    end do

    ! Interior faces
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 > 0) then

        upwd1 = 0.5 * max( v_flux % n(s), 0.0)
        upwd2 = 0.5 * max(-v_flux % n(s), 0.0)

        A % val(A % dia(c1)) = A % val(A % dia(c1)) + upwd1
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + upwd2

        A % val(A % pos(1,s)) =  - upwd2
        A % val(A % pos(2,s)) =  - upwd1

        b(c1) = b(c1) - ( upwd1 * fun % o(c1) - upwd2 * fun % o(c2) )
        b(c2) = b(c2) - ( upwd2 * fun % o(c2) - upwd1 * fun % o(c1) )

      end if
    end do

  else if(fun % adv_scheme .eq. CICSAM .or. &
          fun % adv_scheme .eq. STACS) then

    ! At boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then

        if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) then
          A % val(A % dia(c1)) = A % val(A % dia(c1)) + v_flux % n(s)
        else
          b(c1) = b(c1) - v_flux % n(s) * fun % n(c2)
        end if

      end if
    end do

    ! Interior faces
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
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

        b(c1) = b(c1) - (upwd1 + upwd3) * fun % o(c1) + upwd1 * fun % o(c2)
        b(c2) = b(c2) - (upwd2 - upwd3) * fun % o(c2) + upwd2 * fun % o(c1)

      end if
    end do

  end if

  !------------------------------------------------!
  !   Calculate main coefficient and source term   !
  !        (Part 2 : temporal contribution)        !
  !------------------------------------------------!

  ! Two time levels; linear interpolation
  if(fun % td_scheme .eq. LINEAR) then
    do c = 1, Grid % n_cells
      a0 = Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + a0
      b(c) = b(c) + a0 * fun % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(fun % td_scheme .eq. PARABOLIC) then
    do c = 1, Grid % n_cells
      a0 = Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0 * a0 * fun % o(c) - 0.5 * a0 * fun % oo(c)
    end do
  end if

  !--------------------------------!
  !   Source due to phase change   !
  !--------------------------------!
  call Vof % Mass_Transfer_Estimate()
  call Vof % Mass_Transfer_Vof_Source(b)

  end subroutine
