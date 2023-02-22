!==============================================================================!
  subroutine Inertial_Term_Momentum(Flow, phi, A, b, dt)
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
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Var_Type)            :: phi
  type(Matrix_Type)         :: A
  real                      :: b(:)
  real                      :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  real                      :: a0
  integer                   :: c, c1, c2, s
  real                      :: xc1, yc1, zc1, xc2, yc2, zc2, dotprod
  real, contiguous, pointer :: dens_f(:), dens_cx(:), dens_cy(:), dens_cz(:)
!==============================================================================!

  call Work % Connect_Real_Face(dens_f)
  call Work % Connect_Real_Cell(dens_cx, dens_cy, dens_cz)

  ! Take alias to Grid
  Grid => Flow % pnt_grid

  ! Compute densities at faces (the same way as in Field_Mod_Buoyancy_Forces)
  call Field_Mod_Interpolate_To_Faces_Harmonic(Flow, dens_f, Flow % density)

    ! Calculate densities at cells
    ! Very similar to what you are doing in Field_Mod_Buoyancy_Forces
    dens_cx(:) = 0.0
    dens_cy(:) = 0.0
    dens_cz(:) = 0.0

    do s = 1, Grid % n_faces
      c1  = Grid % faces_c(1, s)
      c2  = Grid % faces_c(2, s)
      xc1 = Grid % xc(c1)
      yc1 = Grid % yc(c1)
      zc1 = Grid % zc(c1)

      ! Units in next three expressions which follow
      ! kg/m^3 * m^2 * m = kg
      dens_cx(c1) = dens_cx(c1)          &
                  + dens_f(s)            &
                  * Grid % sx(s) * (Grid % xf(s) - xc1)

      dens_cy(c1) = dens_cy(c1)          &
                  + dens_f(s)            &
                  * Grid % sy(s) * (Grid % yf(s) - yc1)

      dens_cz(c1) = dens_cz(c1)          &
                  + dens_f(s)            &
                  * Grid % sz(s) * (Grid % zf(s) - zc1)

      if(c2 > 0) then
        xc2 = Grid % xc(c1) + Grid % dx(s)
        yc2 = Grid % yc(c1) + Grid % dy(s)
        zc2 = Grid % zc(c1) + Grid % dz(s)

        ! Units in next three expressions which follow
        ! kg/m^3 * m^2 * m = kg
        dens_cx(c2) = dens_cx(c2)          &
                    - dens_f(s)            &
                    * Grid % sx(s) * (Grid % xf(s) - xc2)

        dens_cy(c2) = dens_cy(c2)          &
                    - dens_f(s)            &
                    * Grid % sy(s) * (Grid % yf(s) - yc2)

        dens_cz(c2) = dens_cz(c2)          &
                    - dens_f(s)            &
                    * Grid % sz(s) * (Grid % zf(s) - zc2)
      end if
    end do

    ! Correct units for density
    do c = 1, Grid % n_cells
      dens_cx(c) = dens_cx(c) / Grid % vol(c)
      dens_cy(c) = dens_cy(c) / Grid % vol(c)
      dens_cz(c) = dens_cz(c) / Grid % vol(c)
    end do


  ! Two time levels; Linear interpolation
  if(phi % td_scheme .eq. LINEAR) then
    do c = 1, Grid % n_cells
      a0 = Flow % density(c) * Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(phi % td_scheme .eq. PARABOLIC) then
    do c = 1, Grid % n_cells
      a0 = Flow % density(c) * Grid % vol(c) / dt
      A % val(A % dia(c)) = A % val(A % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  call Work % Disconnect_Real_Face(dens_f)
  call Work % Disconnect_Real_Cell(dens_cx, dens_cy, dens_cz)

  end subroutine
