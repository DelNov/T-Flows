!==============================================================================!
  subroutine Face_Mod_Allocate(phi, Grid, name_phi)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to allocate memory for a face-centered variable
!>  in a computational grid. It sets up the variable phi with new (n), old (o),
!>  and older (oo) values, as well as an average guessed value (avg). Each of
!>  these arrays is allocated to the size of the total number of faces in the
!>  given grid, and their values are initially set to zero.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Face_Type)         :: phi       !! parent Face_Type
  type(Grid_Type), target :: Grid      !! grid on the variable is defined
  character(len=*)        :: name_phi  !! name of the variable
!==============================================================================!

  ! Store grid for which the variable is defined
  phi % pnt_grid => Grid

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) and old (o) time step
  allocate (phi % n (Grid % n_faces));  phi % n  = 0.0
  allocate (phi % o (Grid % n_faces));  phi % o  = 0.0
  allocate (phi % oo(Grid % n_faces));  phi % oo = 0.0

  ! Average guessed value, guessed value
  allocate (phi % avg (Grid % n_faces));  phi % avg = 0.0

  end subroutine
