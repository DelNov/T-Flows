!==============================================================================!
  subroutine Calculate_Mass_Field_Rate(flow, mult)
!------------------------------------------------------------------------------!
!   Calculate mass flow rate at cell faces based on velocities only.           !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,       only: Grid_Type
  use Field_Mod,      only: Field_Type, Field_Mod_Alias_Momentum, density
  use Var_Mod,        only: Var_Type
  use Multiphase_Mod, only: Multiphase_Type, phase_dens, multiphase_model,  &
                            VOLUME_OF_FLUID
  use Work_Mod,       only: dens_face => r_face_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  real,            pointer :: flux(:)
  integer                  :: s, c1, c2
  real                     :: fs
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  if (multiphase_model .eq. VOLUME_OF_FLUID) then
    do s=1, grid % n_faces
      dens_face(s) = mult % vof_f(s)         * phase_dens(1)  &
                   + (1.0 - mult % vof_f(s)) * phase_dens(2)
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
