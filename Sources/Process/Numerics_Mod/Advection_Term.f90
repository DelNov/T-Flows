!==============================================================================!
  subroutine Numerics_Mod_Advection_Term(phi, coef, v_flux, b)
!------------------------------------------------------------------------------!
!   Purpose: Dicretize advection term in conservation equations.               !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Var_Type) :: phi
  real           :: coef(-phi % pnt_grid % n_bnd_cells:  &
                          phi % pnt_grid % n_cells)
  real           :: v_flux(phi % pnt_grid % n_faces)
  real           :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  real                      :: phif          ! phi and coef at the cell face
  integer                   :: c, c1, c2, s
  real, contiguous, pointer :: phi_min(:), phi_max(:), advect(:), upwind(:)
!==============================================================================!

  ! Take alias to Grid
  Grid => phi % pnt_grid

  call Work % Connect_Real_Cell(phi_min, phi_max, advect, upwind)

  !----------------------------------------------------------------------------!
  !   Compute phi % max and phi % min (used in Numerics_Mod_Advection_Scheme)  !
  !----------------------------------------------------------------------------!
  if(phi % adv_scheme .ne. CENTRAL) then
    call Numerics_Mod_Min_Max(phi, phi_min, phi_max)
  end if

  !----------------!
  !   New values   !
  !----------------!
  advect(:) = 0.0
  upwind(:) = 0.0

  !----------------------------------!
  !   Browse through all the faces   !
  !----------------------------------!
  do s=1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! This could be computed with gradient extrapolation
    phif =      Grid % f(s)  * phi % n(c1)   &
         + (1.0-Grid % f(s)) * phi % n(c2)

    ! Compute phif with desired advection scheme
    if(phi % adv_scheme .ne. CENTRAL) then
      call Numerics_Mod_Advection_Scheme(phif, s, phi, phi_min, phi_max, v_flux)
    end if

    ! Compute advection term (volume-conservative form)
    if(c2 > 0) then
      advect(c1) = advect(c1) - v_flux(s) * phif * coef(c1)
      advect(c2) = advect(c2) + v_flux(s) * phif * coef(c2)
    else
      advect(c1) = advect(c1) - v_flux(s) * phif * coef(c1)
    end if

    ! Store upwinded part of the advection term
    if(v_flux(s) .lt. 0) then   ! from c2 to c1
      upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c2) * coef(c1)
      if(c2 > 0) then
        upwind(c2) = upwind(c2) + v_flux(s) * phi % n(c2) * coef(c2)
      end if
    else
      upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c1) * coef(c1)
      if(c2 > 0) then
        upwind(c2) = upwind(c2) + v_flux(s) * phi % n(c1) * coef(c2)
      end if
    end if

  end do  ! through faces

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, Grid % n_cells
    b(c) = b(c) + advect(c) - upwind(c)
  end do

  call Work % Disconnect_Real_Cell(phi_min, phi_max, advect, upwind)

  end subroutine
