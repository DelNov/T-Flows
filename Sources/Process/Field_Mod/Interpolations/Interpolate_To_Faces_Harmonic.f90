!==============================================================================!
  subroutine Interpolate_To_Faces_Harmonic(Flow, phi_f, phi_c)
!------------------------------------------------------------------------------!
!   Calculates harmonic average on all faces in the domain, from the values    !
!   in the cell centers which surround it.  Useful for density interpolation   !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phi_f(  Flow % pnt_grid % n_faces)
  real                      :: phi_c( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2, reg
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Refresh buffers for gradient components was here, but it is not needed

  ! Perform harmonic average for boundary faces
  ! (Why doesn't it take care of boundary conditions? - check this!)
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      phi_f(s) = phi_c(c1)
    end do
  end do

  ! Perform harmonic average for inside faces
  do s = Faces_In_Domain_And_At_Buffers()
    c1  = Grid % faces_c(1,s)
    c2  = Grid % faces_c(2,s)

    phi_f(s) = Math % Harmonic_Mean(phi_c(c1), phi_c(c2))
  end do

  end subroutine
