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
  type(Grid_Type), pointer :: grid
  real                     :: phi_f          ! phi and coef at the cell face
  integer                  :: c, c1, c2, s
!==============================================================================!

  ! Take alias to grid
  grid => phi % pnt_grid

  !----------------------------------------------------------------------------!
  !   Compute phi % max and phi % min (used in Numerics_Mod_Advection_Scheme)  !
  !----------------------------------------------------------------------------!
  if(phi % adv_scheme .ne. CENTRAL) then
    call Numerics_Mod_Min_Max(phi)
  end if

  !----------------!
  !   New values   !
  !----------------!
  do c = 1, grid % n_cells
    phi % a(c) = 0.0
    phi % c(c) = 0.0  ! use phi % c for upwind advective fluxes
  end do

  !----------------------------------!
  !   Browse through all the faces   !
  !----------------------------------!
  do s=1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! This could be computed with gradient extrapolation
    phi_f =      grid % f(s)  * phi % n(c1)   &
          + (1.0-grid % f(s)) * phi % n(c2)

    ! Compute phi_f with desired advection scheme
    if(phi % adv_scheme .ne. CENTRAL) then
      call Numerics_Mod_Advection_Scheme(phi_f, s, phi, v_flux)
    end if

    ! Compute advection term (non-conservative form)
    if(c2 > 0) then
      phi % a(c1) = phi % a(c1) - v_flux(s) * phi_f * coef(c1)
      phi % a(c2) = phi % a(c2) + v_flux(s) * phi_f * coef(c2)
    else
      phi % a(c1) = phi % a(c1) - v_flux(s) * phi_f * coef(c1)
    end if

    ! Store upwinded part of the advection term in "c"
    if(v_flux(s) .lt. 0) then   ! from c2 to c1
      phi % c(c1) = phi % c(c1) - v_flux(s) * phi % n(c2) * coef(c1)
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2) + v_flux(s) * phi % n(c2) * coef(c2)
      end if
    else
      phi % c(c1) = phi % c(c1) - v_flux(s) * phi % n(c1) * coef(c1)
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2) + v_flux(s) * phi % n(c1) * coef(c2)
      end if
    end if

  end do  ! through faces

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + phi % a(c) - phi % c(c)
    phi % c(c) = 0.0                       ! set it back to zero
  end do

  end subroutine
