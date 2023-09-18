!==============================================================================!
  subroutine Numerics_Mod_Advection_Term(phi, coef, v_flux, b, blend_matrix)
!------------------------------------------------------------------------------!
!   Purpose: Dicretize advection term in conservation equations.               !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Var_Type), intent(in)     :: phi
  real,           intent(in)     :: coef(-phi % pnt_grid % n_bnd_cells:  &
                                          phi % pnt_grid % n_cells)
  real,           intent(in)     :: v_flux(phi % pnt_grid % n_faces)
  real,           intent(inout)  :: b(:)
  logical,        intent(in)     :: blend_matrix
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer  :: Grid
  real                      :: phif          ! phi and coef at the cell face
  integer                   :: c, c1, c2, s, reg
  real, contiguous, pointer :: phi_min(:), phi_max(:), advect(:), upwind(:)
!==============================================================================!

  call Profiler % Start('Advection_Term')

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

  !----------------------------------------------!
  !   Compute advection term on the boundaries   !
  !----------------------------------------------!
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! This could be computed with gradient extrapolation
      phif =      Grid % f(s)  * phi % n(c1)   &
           + (1.0-Grid % f(s)) * phi % n(c2)

      ! Compute phif with desired advection scheme
      if(phi % adv_scheme .ne. CENTRAL) then
        call Numerics_Mod_Advection_Scheme(phif, s, phi,  &
                                           phi_min, phi_max, v_flux)
      end if

      ! Compute advection term (volume-conservative form)
      advect(c1) = advect(c1) - v_flux(s) * phif * coef(c1)

    end do  ! faces
  end do    ! regions

  !-------------------------------------------!
  !   Compute upwind term on the boundaries   !
  !-------------------------------------------!
  if(blend_matrix) then
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! Store upwinded part of the advection term
        if(v_flux(s) .lt. 0) then   ! from c2 to c1
          upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c2) * coef(c1)
        else
          upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c1) * coef(c1)
        end if

      end do  ! faces
    end do    ! regions
  end if

  !----------------------------------------------!
  !   Compute advection term inside the domain   !
  !----------------------------------------------!
  do s = Faces_In_Domain()
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
    advect(c1) = advect(c1) - v_flux(s) * phif * coef(c1)
    advect(c2) = advect(c2) + v_flux(s) * phif * coef(c2)

  end do  ! through faces in domain

  !-------------------------------------------!
  !   Compute upwind term inside the domain   !
  !-------------------------------------------!
  if(blend_matrix) then
    do s = Faces_In_Domain()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! Store upwinded part of the advection term
      if(v_flux(s) .lt. 0) then   ! from c2 to c1
        upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c2) * coef(c1)
        upwind(c2) = upwind(c2) + v_flux(s) * phi % n(c2) * coef(c2)
      else
        upwind(c1) = upwind(c1) - v_flux(s) * phi % n(c1) * coef(c1)
        upwind(c2) = upwind(c2) + v_flux(s) * phi % n(c1) * coef(c2)
      end if

    end do  ! through faces in domain
  end if

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = Cells_In_Domain()
    b(c) = b(c) + advect(c) - upwind(c)
  end do

  call Work % Disconnect_Real_Cell(phi_min, phi_max, advect, upwind)

  call Profiler % Stop('Advection_Term')

  end subroutine
