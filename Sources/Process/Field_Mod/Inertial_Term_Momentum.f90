!==============================================================================!
  subroutine Field_Mod_Inertial_Term_Momentum(flow, phi, a, b, dt)
!------------------------------------------------------------------------------!
!   Purpose: Dicretize inertial term in momentum conservation equations.       !
!                                                                              !
!   To my utter surprise: this doesn't seem to be needed for force balance     !
!   on the faces.  I don't understand why - I wanted to compute inertial       !
!   terms on the faces with densities which are the same as forces on the      !
!   faces in Rhie and Chow - but it doesn't seem to be needed.  Luckily so,    !
!   I dare to say, because if this was needed we would have three different    !
!   systems of linear equations for each velocity component - no idea how      !
!   would I use them to form pressure correction matrices :-o                  !
!                                                                              !
!   Yet, it is surprising and I better keep it for a while longer in the       !
!   code for a subsequent analysis - should I ever have time for that.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: dens_f  => r_face_01,  &
                      dens_cx => r_cell_01,  &
                      dens_cy => r_cell_02,  &
                      dens_cz => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type)           :: phi
  type(Matrix_Type)        :: a
  real                     :: b(:)
  real                     :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: a0
  integer                  :: c, c1, c2, s
  real                     :: xc1, yc1, zc1, xc2, yc2, zc2, dotprod
!==============================================================================!

  ! Take alias to grid
  grid => flow % pnt_grid

  ! Compute densities at faces (the same way as in Field_Mod_Buoyancy_Forces)
  call Field_Mod_Interpolate_To_Faces_Harmonic(flow, dens_f, flow % density)

    ! Calculate densities at cells
    ! Very similar to what you are doing in Field_Mod_Buoyancy_Forces
    dens_cx(:) = 0.0
    dens_cy(:) = 0.0
    dens_cz(:) = 0.0

    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)
      xc1 = grid % xc(c1)
      yc1 = grid % yc(c1)
      zc1 = grid % zc(c1)

      ! Units in next three expressions which follow
      ! kg/m^3 * m^2 * m = kg
      dens_cx(c1) = dens_cx(c1)          &
                  + dens_f(s)            &
                  * grid % sx(s) * (grid % xf(s) - xc1)

      dens_cy(c1) = dens_cy(c1)          &
                  + dens_f(s)            &
                  * grid % sy(s) * (grid % yf(s) - yc1)

      dens_cz(c1) = dens_cz(c1)          &
                  + dens_f(s)            &
                  * grid % sz(s) * (grid % zf(s) - zc1)

      if(c2 > 0) then
        xc2 = grid % xc(c1) + grid % dx(s)
        yc2 = grid % yc(c1) + grid % dy(s)
        zc2 = grid % zc(c1) + grid % dz(s)

        ! Units in next three expressions which follow
        ! kg/m^3 * m^2 * m = kg
        dens_cx(c2) = dens_cx(c2)          &
                    - dens_f(s)            &
                    * grid % sx(s) * (grid % xf(s) - xc2)

        dens_cy(c2) = dens_cy(c2)          &
                    - dens_f(s)            &
                    * grid % sy(s) * (grid % yf(s) - yc2)

        dens_cz(c2) = dens_cz(c2)          &
                    - dens_f(s)            &
                    * grid % sz(s) * (grid % zf(s) - zc2)
      end if
    end do

    ! Correct units for density
    do c = 1, grid % n_cells
      dens_cx(c) = dens_cx(c) / grid % vol(c)
      dens_cy(c) = dens_cy(c) / grid % vol(c)
      dens_cz(c) = dens_cz(c) / grid % vol(c)
    end do


  ! Two time levels; Linear interpolation
  if(phi % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = flow % density(c) * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(phi % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = flow % density(c) * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  end subroutine
