!==============================================================================!
  subroutine Field_Mod_Interpolate_To_Faces_Harmonic(flow, phi_f, phi_c)
!------------------------------------------------------------------------------!
!   Calculates harmonic average on all faces in the domain, from the values    !
!   in the cell centers which surround it.  Useful for density interpolation   !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phi_f(  flow % pnt_grid % n_faces)
  real                     :: phi_c( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers for gradient components was here, but it is not needed

  ! Perform harmonic average on all faces
  do s = 1, grid % n_faces
    c1  = grid % faces_c(1, s)
    c2  = grid % faces_c(2, s)

    if(c2 > 0) then
      phi_f(s) = Math_Mod_Harmonic_Mean(phi_c(c1), phi_c(c2))
    else
      phi_f(s) = phi_c(c1)
    end if

    phi_f(s) = phi_f(s)
  end do

  end subroutine
