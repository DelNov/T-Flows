!==============================================================================!
  subroutine Turb_Mod_Src_F22_Rsm_Manceau_Hanjalic(turb, sol)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation equation                      !
!   and imposing  boundary condition forf22                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: f22
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: sor11,  f22hg
!==============================================================================!
!                                                                              !
!   The form of source terms are :                                             !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22hg*dV ; f22hg - f22hg homogenious is placed in a source              !
!    |                     coefficients b(c)                                   !
!   /                                                                          !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/kin(c) - 2.0/3.0)/Tsc(c)     &             !
!              + 2.0*Cv_2*p_kin(c)/(3.0*kin(c))                                !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22*dV ;this term is placed in a diagonal of coefficient matrice        !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!                                                                              !
!  Dimensions of certain variables                                             !
!                                                                              !
!     Tsc            [s]                                                       !
!     kin            [m^2/s^2]                                                 !
!     eps            [m^3/s^2]                                                 !
!     vi2            [m^2/s^2]                                                 !
!     f22            [-]                                                       !
!     l_scale        [m]                                                       !
!                                                                              !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  f22  => turb % f22
  call Solver_Mod_Alias_System(sol, a, b)

  call Time_And_Length_Scale(grid, turb)

  ! Source term f22hg
  do c = 1, grid % n_cells
    f22hg = 1.0
    sor11 = grid % vol(c) / turb % l_scale(c)**2
    a % val(a % dia(c)) = a % val(a % dia(c)) + sor11
    b(c) = b(c) + f22hg * grid % vol(c) / turb % l_scale(c)**2
  end do

  ! Source term
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

          f22 % n(c2) = 0.0

      end if   ! end if of BC=wall
    end if    ! end if of c2<0
  end do

  end subroutine
