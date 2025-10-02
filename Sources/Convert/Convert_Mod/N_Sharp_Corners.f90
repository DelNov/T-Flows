!==============================================================================!
  integer function N_Sharp_Corners(Convert, Grid, sharp_corner)
!------------------------------------------------------------------------------!
!>  This function is designed to identify and count sharp corners in a 3D
!>  grid structure. It takes the Grid as an input and returns an array
!>  sharp_corner marking the nodes at sharp corners.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: The function initializes a counter cnt and sets all      !
!     entries in the sharp_corner array to zero.                               !
!   * Identifying Sharp Corners:                                               !
!     - Iterates over each edge in the grid.                                   !
!     - For edges between two boundary faces, it computes the normals of these !
!       faces.                                                                 !
!     - If the angle between these normals is less than 45 degrees (indicating !
!       a sharp corner), the nodes at either end of the edge are marked in     !
!       sharp_corner by incrementing their respective counts.                  !
!   * Counting Sharp Corners:                                                  !
!     - After marking, it counts the number of nodes that have been marked     !
!       more than twice, indicating a sharp corner.                            !
!     - The count is updated in the sharp_corner array, while nodes not part   !
!       of a sharp corner are reset to zero.                                   !
!   * Return Value:                                                            !
!     - The function returns the total count of sharp corner nodes.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  type(Grid_Type)     :: Grid
  integer             :: sharp_corner(Grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: e, cnt, s1, s2, n1, n2, n
  real    :: norm_1(3), norm_2(3)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  cnt             = 0
  sharp_corner(:) = 0

  !---------------------------------------------!
  !   Fetch geometrically sharp edges (first)   !
  !---------------------------------------------!
  do e = 1, Grid % n_edges

    ! If edge is in-between two boundary faces
    ! (yes, Grid % edges_f stores only boundary faces, very bad name choice)
    if(Grid % edges_fb(1, e) .ne. 0 .and.  &
       Grid % edges_fb(2, e) .ne. 0) then

      ! Take the boundary faces around the edge
      s1 = Grid % edges_fb(1, e)
      s2 = Grid % edges_fb(2, e)

      ! Compute normals of the two faces surrounding the edge
      norm_1(1) = Grid % sx(s1);  norm_2(1) = Grid % sx(s2)
      norm_1(2) = Grid % sy(s1);  norm_2(2) = Grid % sy(s2)
      norm_1(3) = Grid % sz(s1);  norm_2(3) = Grid % sz(s2)
      norm_1(1:3) = norm_1(1:3) / norm2(norm_1(1:3))
      norm_2(1:3) = norm_2(1:3) / norm2(norm_2(1:3))

      ! If angle is less than 45 (135) degrees, or is simply sharp
      if(dot_product(norm_1(1:3), norm_2(1:3)) < 0.7071) then

        ! Take edges' nodes and mark them (increase
        ! the number of times they have been visited)
        n1 = Grid % edges_n(1, e)
        n2 = Grid % edges_n(2, e)

        sharp_corner(n1) = sharp_corner(n1) + 1
        sharp_corner(n2) = sharp_corner(n2) + 1
      end if
    end if

  end do

  !------------------------------------------------------------!
  !   Count nodes which have been marked more than two times   !
  !------------------------------------------------------------!
  do n = 1, Grid % n_nodes
    if(sharp_corner(n) > 2) then
      cnt = cnt + 1
      sharp_corner(n) = cnt
    else
      sharp_corner(n) = 0
    end if
  end do

  ! Return a sum of all marked nodes
  N_Sharp_Corners = cnt

  end function

