!==============================================================================!
  subroutine Source_F22_K_Eps_Zeta_F(flow, sol)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation  equation                     !
!   for vi2 and imposing  boundary condition for f22                           !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Grid_Mod,   only: Grid_Type
  use Field_Mod
  use Turb_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: sor_11, f22hg
!==============================================================================!
!                                                                              !
!  The form of source terms are :                                              !
!                                                                              !
!      int( f22hg*dV ),                                                        !
!  where f22hg - f22hg homogenious is placed in a source coefficients b(c)     !
!                                                                              !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/kin(c) - 2.0/3.0)/t_scale(c)     &         !
!              + 2.0*Cv2*p_kin(c)/(3.0*kin(c))                                 !
!                                                                              !
!    int( f22*dV ); this term is placed in a diagonal of coefficient matrix    !
!                                                                              !
!                                                                              !
!  Dimensions of certain variables                                             !
!                                                                              !
!     t_scale        [s]                                                       !
!     kin            [m^2/s^2]                                                 !
!     eps            [m^3/s^2]                                                 !
!     vi2            [m^2/s^2]                                                 !
!     f22            [-]                                                       !
!     l_scale        [m]                                                       !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  a    => sol % a
  b    => sol % b % val

  call Time_And_Length_Scale(grid)

 ! Source term f22hg
 do c = 1, grid % n_cells
   f22hg = (1.0 - c_f1 - 0.65 * p_kin(c)/density  &
         / (eps  % n(c) + TINY))                  &
         * (zeta % n(c) - TWO_THIRDS)             &
         / (t_scale(c) + TINY)                    &
         + 0.0085 * (p_kin(c)/density) / (kin % n(c) + TINY)
   b(c) = b(c) + f22hg * grid % vol(c) / (l_scale(c)**2 + TINY)
 end do

 ! Source term f22hg
 do c = 1, grid % n_cells
   sor_11 = grid % vol(c)/(l_scale(c)**2 + TINY)
   a % val(a % dia(c)) = a % val(a % dia(c)) + sor_11 
 end do

 ! Imposing boundary condition for f22 on the wall
 do s = 1, grid % n_faces
   c1 = grid % faces_c(1,s)
   c2 = grid % faces_c(2,s)
   if(c2 < 0) then
     if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
        Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

       f22 % n(c2) = -2.0 * viscosity/density * zeta % n(c1) &
                   / grid % wall_dist(c1)**2

      ! Fill in a source coefficients

      ! Linearization of the near wall terms helps to get more  
       ! stable solution, especially for small wall distance.
       b(c1) = b(c1) + a % fc(s) * f22 % n(c2)
     end if   ! end if of BC=wall
   end if    ! end if of c2<0
 end do

 end subroutine
