!==============================================================================!
  subroutine Source_Eps_K_Eps(grid)
!------------------------------------------------------------------------------!
!   Computes the source terms in the eps transport equation,                   !
!   wall shear stress (wall function approuch)                                 !
!------------------------------------------------------------------------------!
!   int( density (c_1e eps/kin Gk - c_2e eps^2/kin) )dV                        !
!                                                                              !
!   assigns epsilon from the wall function:                                    !
!                                                                              !
!   Eps_w = Cmu^(3/4)* Kin^(3/2)/(Kappa/y)                                     !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
  use Work_Mod, only: shear_x => r_cell_01,  &
                      shear_y => r_cell_02,  &
                      shear_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, j
  real              :: re_t, f_mu, l1, l2, yap
!==============================================================================!
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

  if(turbulence_model_variant .eq. HIGH_RE) then
    do c = 1, grid % n_cells
      ! Positive contribution:
      b(c) = b(c) + &
        c_1e * density * eps % n(c) / kin % n(c) * p_kin(c) * grid % vol(c)

      ! Negative contribution:
      A % val(A % dia(c)) = A % val(A % dia(c)) &
           + c_2e * density * eps % n(c) / kin % n(c) * grid % vol(c)
    end do

    !--------------------------------------------!
    !   Cut-off the wall influence in order to   !
    !   impose the boundary condition for EPS    !
    !--------------------------------------------!
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

          ! This will fix the value of eps in the first cell
          if(rough_walls .eq. NO) then
            eps % n(c1) = c_mu75 * (kin % n(c1))**1.5   &
                        / (kappa*grid % wall_dist(c1))
          else if(rough_walls .eq. YES) then
            eps % n(c1) = c_mu75*(kin % n(c1))**1.5  &
                        / (kappa*(grid % wall_dist(c1)+Zo))
          end if

          do j = A % row(c1), A % row(c1+1) -1
            A % val(j) = 0.0
          end do

          A % val(A % dia(c1)) = 1.0 * density
          b(c1) = density * eps % n(c1) * grid % vol(c1)

        end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
      end if    ! c2 < 0
    end do
  end if        ! end if mode = wf

  if(turbulence_model_variant .eq. LOW_RE) then
  !-------------------------!
  !   Jones-Launder model   !
  !-------------------------!

   call Calculate_Shear_And_Vorticity(grid)
   call Grad_Mod_For_Phi(grid, shear, 1, shear_x, .true.)  ! dU/dx
   call Grad_Mod_For_Phi(grid, shear, 2, shear_y, .true.)  ! dW/dy
   call Grad_Mod_For_Phi(grid, shear, 3, shear_z, .true.)  ! dV/dz

    do c = 1, grid % n_cells
      ! Positive contribution:
      b(c) = b(c) + &
        c_1e * density * eps % n(c) / kin % n(c) * p_kin(c) * grid % vol(c) &
        + 2.0 * viscosity * vis_t(c) / density * &
           (shear_x(c)**2. + shear_y(c)**2. + shear_z(c)**2.) * grid % vol(c)

      ! Negative contribution:
      re_t = kin % n(c)**2./((viscosity/density)*eps % n(c))
      f_mu = 1.0 - 0.3*exp(-(re_t**2.))
      A % val(A % dia(c)) = A % val(A % dia(c)) &
        + f_mu * c_2e * density * eps % n(c) / kin % n(c) * grid % vol(c)

       ! Yap correction
       l1 = kin % n(c)**1.5/eps % n(c)
       l2 = 2.55 * grid % wall_dist(c)
       yap = 0.83 * eps % n(c)**2./kin % n(c)  &
                  * max((l1/l2 - 1.0) * (l1/l2)**2., 0.0)
       b(c) = b(c) + yap * density* grid % vol(c)
    end do

    !-----------------------------------!
    !   Boundary condition fo epsilon   !
    !-----------------------------------!
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

          eps % n(c2) = 0.0

        end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
      end if    ! c2 < 0
    end do
  end if

  end subroutine
