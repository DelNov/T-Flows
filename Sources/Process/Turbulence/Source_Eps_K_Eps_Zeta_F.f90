!==============================================================================!
  subroutine Source_Eps_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Calculates source terms in equation of dissipation of turbulent energy     !
!   and imposes boundary condition                                             !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  integer :: c, s, c1, c2,j 
  real    :: Esor, c_11e, EBF
  real    :: EpsWall, EpsHom
  real    :: y_plu
!==============================================================================!
!   In dissipation of turbulent kinetic energy equation exist two              !
!   source terms which have form:                                              !
!                                                                              !
!    int( density ((Cv_e1 * p_kin - Cv_11 eps) / t_scale) * dV                 !
!                                                                              !
!   First, positive , source term is solved and added to source  coefficient   !
!   b(c) on right hand side.  Second, negative, source term is added to main   !
!   diagonal left hand side coefficient matrix in order to increase stability  !
!   of solver.  It is nessesary to calculate coefficient Cv_11 using kin,      !
!   Cv_e2, vi2 and coefficient A1                                              !
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!
  call Time_And_Length_Scale(grid)

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    do c = 1, grid % n_cells 
      Esor = grid % vol(c)/(t_scale(c)+TINY)
      c_11e = c_1e*(1.0 + alpha * ( 1.0/(zeta % n(c)+TINY) ))    
      b(c) = b(c) + c_11e * density * Esor * p_kin(c)
 
      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + c_2e * Esor * density
    end do                   
  end if

  !-------------------------------------------------------!
  !   Following block shows density dependent behaviour   !
  !-------------------------------------------------------!

  ! Imposing a boundary condition on wall for eps
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER ) then
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
          Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        EpsWall = 2.0*viscosity/density * kin % n(c1) / grid % wall_dist(c1)**2
        EpsHom = c_mu75 * kin % n(c1)**1.5 / (grid % wall_dist(c1) * kappa)
        u_tau(c1) = c_mu25 * kin % n(c1)**0.5

        y_plu = c_mu25 * sqrt(kin % n(c1)) * grid % wall_dist(c1) / &
          (viscosity/density) ! standard
        EBF = 0.001*y_plu**4.0/(1.0 + y_plu)
        eps % n(c1) = EpsWall * exp(-1.0 * EBF) + EpsHom * exp(-1.0 / EBF)
        
        if(rough_walls .eq. YES) then
          eps % n(c1) = c_mu75*kin % n(c1)**1.5 / (grid % wall_dist(c1) * kappa)
        end if

        ! Adjusting coefficient to fix eps value in near wall calls
        do j = A % row(c1), A % row(c1 + 1) - 1
          A % val(j) = 0.0
        end do

        b(c1) = eps % n(c1) * density
        A % val(A % dia(c1)) = 1.0 * density
      end if
    end if
  end do  

  end subroutine
