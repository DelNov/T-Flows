!==============================================================================!
  subroutine Rhie_And_Chow(flow, mult, sol, u_f, v_f, w_f)
!------------------------------------------------------------------------------!
!   Computes face velocitites with Rhie and Chow interpolation method          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  real                          :: u_f(flow % pnt_grid % n_faces)
  real                          :: v_f(flow % pnt_grid % n_faces)
  real                          :: w_f(flow % pnt_grid % n_faces)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c1, c2
!==============================================================================!

  call Cpu_Timer_Mod_Start('Rhie_And_Chow')

  ! Take aliases
  grid    => flow % pnt_grid
  a       => sol % a
  b       => sol % b % val
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  ! call User_Mod_Beginning_Of_Compute_Pressure(flow, mult, ini)

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Face is inside the domain
    if(c2 > 0) then

      ! Interpolate velocity

      ! If there is a jump in velocities, call specialized gradient calculation
      if(flow % mass_transfer) then
        u_f(s) = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, u, s)
        v_f(s) = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, v, s)
        w_f(s) = Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump(mult, w, s)

      ! No jumps, call usual routines
      else
        u_f(s) = Field_Mod_Interpolate_Var_To_Face(flow, u, s)
        v_f(s) = Field_Mod_Interpolate_Var_To_Face(flow, v, s)
        w_f(s) = Field_Mod_Interpolate_Var_To_Face(flow, w, s)
      end if

    end if

  end do

  call Cpu_Timer_Mod_Stop('Rhie_And_Chow')

  end subroutine
