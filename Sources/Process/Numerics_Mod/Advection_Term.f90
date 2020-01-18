!==============================================================================!
  subroutine Numerics_Mod_Advection_Term(phi, coef, flux, sol,  &
                                         phi_i, phi_j, phi_k,   &
                                         di, dj, dk)
!------------------------------------------------------------------------------!
!   Purpose: Dicretize advection term in conservation equations.               !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Var_Type)            :: phi
  real                      :: phi_i(-phi % pnt_grid % n_bnd_cells:  &
                                      phi % pnt_grid % n_cells),     &
                               phi_j(-phi % pnt_grid % n_bnd_cells:  &
                                      phi % pnt_grid % n_cells),     &
                               phi_k(-phi % pnt_grid % n_bnd_cells:  &
                                      phi % pnt_grid % n_cells)
  real                      :: di(phi % pnt_grid % n_faces),         &
                               dj(phi % pnt_grid % n_faces),         &
                               dk(phi % pnt_grid % n_faces)
  real                      :: coef(-phi % pnt_grid % n_bnd_cells:  &
                                     phi % pnt_grid % n_cells)
  real                      :: flux(phi % pnt_grid % n_faces)
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real,            pointer :: b(:)
  real                     :: phi_f, coef_f  ! phi and coef at the cell face
  integer                  :: c, c1, c2, s
!==============================================================================!

  ! Take alias to grid
  grid => phi % pnt_grid
  b    => sol % b % val

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

    phi_f =      grid % f(s)  * phi % n(c1)   &
          + (1.0-grid % f(s)) * phi % n(c2)

    coef_f =      grid % f(s)  * coef(c1)   &
           + (1.0-grid % f(s)) * coef(c2)

    ! Compute phi_f with desired advection scheme
    if(phi % adv_scheme .ne. CENTRAL) then
      call Numerics_Mod_Advection_Scheme(phi_f,    &
                                         s,        &
                                         phi,      &
                                         phi_i,    &
                                         phi_j,    &
                                         phi_k,    &
                                         di,       &
                                         dj,       &
                                         dk,       &
                                         flux)
    end if

    ! Compute advection term
    if(c2 > 0) then
      phi % a(c1) = phi % a(c1) - flux(s) * phi_f * coef_f
      phi % a(c2) = phi % a(c2) + flux(s) * phi_f * coef_f
    else
      phi % a(c1) = phi % a(c1) - flux(s) * phi_f * coef_f
    end if

    ! Store upwinded part of the advection term in "c"
    if(flux(s) .lt. 0) then   ! from c2 to c1
      phi % c(c1) = phi % c(c1) - flux(s) * phi % n(c2) * coef(c2)
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2) + flux(s) * phi % n(c2) * coef(c2)
      end if
    else
      phi % c(c1) = phi % c(c1) - flux(s) * phi % n(c1) * coef(c1)
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2) + flux(s) * phi % n(c1) * coef(c1)
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
