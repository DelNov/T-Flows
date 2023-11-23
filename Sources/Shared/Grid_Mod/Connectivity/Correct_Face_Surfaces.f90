!==============================================================================!
  subroutine Correct_Face_Surfaces(Grid)
!------------------------------------------------------------------------------!
!>  Correct face surface areas and (hopefully) order their nodes.
!------------------------------------------------------------------------------!
!   This is important, and I think also called only from, Convert, in cases    !
!   of polyhedral grid conversion.  When creating dual grids, Convert only     !
!   makes sure that faces have continuous nodes, but is not aware of cells     !
!   c1 and c2, their positions and their future numbers.  That's why it is     !
!   probably easier and meaningful to do it in a separate function.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, n
  real    :: dx, dy, dz
!==============================================================================!

  ! If grid is polyhedral, nodes are in arbitrary order, try to fix it
  if(Grid % polyhedral) then

    ! Do the actual calculation
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      if(c2 .gt. c1 .or. c2 .lt. 0) then
        dx = grid % xc(c2) - grid % xc(c1)
        dy = grid % yc(c2) - grid % yc(c1)
        dz = grid % zc(c2) - grid % zc(c1)
        if(Grid % sx(s)*dx + Grid % sy(s)*dy + Grid % sz(s)*dz .lt. 0.0) then
          Grid % sx(s) = - Grid % sx(s)
          Grid % sy(s) = - Grid % sy(s)
          Grid % sz(s) = - Grid % sz(s)
          n = Grid % faces_n_nodes(s)
          call Sort % Reverse_Order_Int(Grid % faces_n(1:n,s))
        end if
      end if

    end do  ! through faces

  end if  ! grid is polyhedral

  print *, '# Face surfaces corrected !'

  end subroutine
