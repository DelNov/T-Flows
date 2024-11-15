!==============================================================================!
  subroutine Src_F22_K_Eps_Zeta_F(Turb, Sol)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation  equation                     !
!   for vi2 and imposing  boundary condition for f22                           !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!----------------------------------[Locals]------------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: kin, eps, f22, zeta
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, reg
  real                       :: sor_11, f22hg
!------------------------------------------------------------------------------!
!                                                                              !
!  The form of source terms are :                                              !
!                                                                              !
!         /                                                                    !
!        |                                                                     !
!        | f22hg * dV                                                          !
!        |                                                                     !
!       /                                                                      !
!                                                                              !
!  where f22hg - f22hg homogenious is placed in a source coefficients b(c)     !
!                                                                              !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/kin(c) - 2.0/3.0)/t_scale(c)     &         !
!              + 2.0*Cv2*p_kin(c)/(3.0*kin(c))                                 !
!                                                                              !
!    int( f22*dV ); this term is placed in a diagonal of coefficient matrix    !
!                                                                              !
!  Dimensions of certain variables                                             !
!                                                                              !
!     t_scale        [s]                                                       !
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
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Sol % Alias_Native       (A, b)

  call Turb % Time_And_Length_Scale(Grid)

  ! Source term f22hg
  do c = Cells_In_Domain()
    f22hg = (1.0 - Turb % c_f1 - 0.65 * Turb % p_kin(c) / Flow % density(c)  &
          / (eps  % n(c) + TINY))                                            &
          * (zeta % n(c) - TWO_THIRDS)                                       &
          / (Turb % t_scale(c) + TINY)                                       &
          + 0.0085 * (Turb % p_kin(c) / Flow % density(c)) / (kin % n(c) + TINY)
    b(c) = b(c) + f22hg * Grid % vol(c) / (Turb % l_scale(c)**2 + TINY)
  end do

  ! Source term f22hg
  do c = Cells_In_Domain()
    sor_11 = Grid % vol(c)/(Turb % l_scale(c)**2 + TINY)
    A % val(A % dia(c)) = A % val(A % dia(c)) + sor_11 
  end do

  ! Imposing boundary condition for f22 on the wall
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)

        f22 % n(c2) = -2.0 * Flow % viscosity(c1)        &
                    / Flow % density(c1) * zeta % n(c1)  &
                    / Grid % wall_dist(c1)**2

        ! Fill in a source coefficients

        ! Linearization of the near wall terms helps to get more
        ! stable solution, especially for small wall distance.
        b(c1) = b(c1) + A % fc(s) * f22 % n(c2)
      end do    ! faces in regions
    end if      ! region is WALL or WALLFL
  end do        ! through regions

  end subroutine
