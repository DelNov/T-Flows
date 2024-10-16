!==============================================================================!
  subroutine Interpolate_Nodes_To_Faces(Flow, var_node, var_face)
!------------------------------------------------------------------------------!
!   Interpolates a variable from nodes to faces                                !
!   It is an experimental function developed as a help to implement Gaussian   !
!   computation of gradients, and hasn't been thoroughly tested.  It is simple !
!   enough (face value is arithmetic average of it's node values) and it did   !
!   pass some basic functionality tests, but didn't prove itself in long and   !
!   rigorous testing campaigns like many other functions in the code. Parallel !
!   version is trivial, provided that nodal values are properly calculated     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type)  :: Flow
  real, intent(in)   :: var_node(1:Flow % pnt_grid % n_nodes)
  real, intent(out)  :: var_face(1:Flow % pnt_grid % n_faces)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, s, i_nod
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Interpolate from all nodes to faces
  do s = 1, Grid % n_faces

    var_face(s) = 0.0

    ! Loop on nodes
    do i_nod = 1, abs(Grid % faces_n_nodes(s))

      n = Grid % faces_n(i_nod, s)
      var_face(s) = var_face(s) + var_node(n)

    end do

    var_face(s) = var_face(s) / real(Grid % faces_n_nodes(s))

  end do

  end subroutine
