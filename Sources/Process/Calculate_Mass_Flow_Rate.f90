!==============================================================================!
  subroutine Calculate_Mass_Field_Rate(flow)
!------------------------------------------------------------------------------!
!   Calculate mass flow rate at cell faces based on velocities only.           !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type, Field_Mod_Alias_Momentum, density,          &
                       multiphase, density_mult
  use Var_Mod,   only: Var_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  real,            pointer :: flux(:)
  integer                  :: s, c1, c2
  real                     :: fs
  real       , allocatable :: dens_face(:)
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  allocate(dens_face(grid % n_faces))

  if (multiphase) then
    do s=1, grid % n_faces
      dens_face(s) = flow % vof_f(s) * density_mult(1)                    &
                   + (1.0 - flow % vof_f(s)) * density_mult(2)
    end do
  else
    do s=1, grid % n_faces
      dens_face(s) = density(grid % faces_c(1,s))
    end do

  end if
  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if( c2 > 0 ) then

      ! Extract the "centred" pressure terms from cell velocities
      flux(s) = dens_face(s) *                                            &
              (   (fs * u % n(c1) + (1.0-fs) * u % n(c2)) * grid % sx(s)  &
                + (fs * v % n(c1) + (1.0-fs) * v % n(c2)) * grid % sy(s)  &
                + (fs * w % n(c1) + (1.0-fs) * w % n(c2)) * grid % sz(s) )

    ! Face is on the boundary
    else
      flux(s) = dens_face(s) *              &
              (   u % n(c1) * grid % sx(s)  &
                + v % n(c1) * grid % sy(s)  &
                + w % n(c1) * grid % sz(s) )
    end if

  end do

  end subroutine
