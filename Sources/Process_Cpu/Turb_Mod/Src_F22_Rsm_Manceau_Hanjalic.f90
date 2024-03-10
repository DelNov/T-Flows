!==============================================================================!
  subroutine Src_F22_Rsm_Manceau_Hanjalic(Turb, Sol)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation equation                      !
!   and imposing  boundary condition forf22                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: f22
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, reg
  real                       :: sor11, f22hg
!------------------------------------------------------------------------------!
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
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  f22  => Turb % f22
  call Sol % Alias_Native(A, b)

  call Turb % Time_And_Length_Scale(Grid)

  ! Source term f22hg
  do c = Cells_In_Domain()
    f22hg = 1.0
    sor11 = Grid % vol(c) / Turb % l_scale(c)**2
    A % val(A % dia(c)) = A % val(A % dia(c)) + sor11
    b(c) = b(c) + f22hg * Grid % vol(c) / Turb % l_scale(c)**2
  end do

  ! Source term
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        f22 % n(c2) = 0.0

      end do  ! faces
    end if    ! boundary condition type
  end do      ! regions

  end subroutine
