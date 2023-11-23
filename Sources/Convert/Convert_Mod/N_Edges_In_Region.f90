!==============================================================================!
  integer function N_Edges_In_Region(Convert, Grid, bc, edge_data)
!------------------------------------------------------------------------------!
!   Counts and marks (with edge_data) edges in the given boundary color        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  type(Grid_Type)     :: Grid
  integer             :: bc
  integer             :: edge_data(Grid % n_edges)
!-----------------------------------[Locals]-----------------------------------!
  integer :: e, cnt, s1, s2, n1, n2
  real    :: norm_1(3), norm_2(3)
  real    :: xs, ys, zs, xe, ye, ze
  real    :: vec_ef(3)  ! from edge to face mid-point
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  edge_data(:) = 0

  !---------------------------------------------!
  !                                             !
  !   Fetch geometrically sharp edges (first)   !
  !                                             !
  !---------------------------------------------!
  do e = 1, Grid % n_edges

    !------------------------!
    !   Edge is on this BC   !
    !------------------------!
    if(Grid % edges_bc(bc, e) .gt. 0) then
    if(Grid % edges_fb(1, e) .ne. 0 .and.  &
       Grid % edges_fb(2, e) .ne. 0) then

      ! Take face indices
      s1 = Grid % edges_fb(1, e)
      s2 = Grid % edges_fb(2, e)

      ! Compute normals of the two faces surrounding the edge
      norm_1(1) = Grid % sx(s1);  norm_2(1) = Grid % sx(s2)
      norm_1(2) = Grid % sy(s1);  norm_2(2) = Grid % sy(s2)
      norm_1(3) = Grid % sz(s1);  norm_2(3) = Grid % sz(s2)
      norm_1(1:3) = norm_1(1:3) / norm2(norm_1(1:3))
      norm_2(1:3) = norm_2(1:3) / norm2(norm_2(1:3))

      !----------------------------------------!
      !   Edge is not sharp, mark it as zero   !
      !----------------------------------------!
      if(dot_product(norm_1, norm_2) >= 0.7071) then
        edge_data(e) = 0

      !--------------------------------------------------!
      !   Edge is geometrically sharp, still check if    !
      !   is in convex or concave corner of a new face   !
      !--------------------------------------------------!
      else
        cnt = cnt + 1

        ! Calculate edges coordinates
        n1 = Grid % edges_n(1, e)
        n2 = Grid % edges_n(2, e)
        xe = 0.5 * (Grid % xn(n1) + Grid % xn(n2))
        ye = 0.5 * (Grid % yn(n1) + Grid % yn(n2))
        ze = 0.5 * (Grid % zn(n1) + Grid % zn(n2))

        ! Calculate midpoint between cell faces surrounding the edge
        xs = 0.5 * (Grid % xf(s1) + Grid % xf(s2))
        ys = 0.5 * (Grid % yf(s1) + Grid % yf(s2))
        zs = 0.5 * (Grid % zf(s1) + Grid % zf(s2))

        ! Vector connecting edge center with face mid-point
        vec_ef(1) = xs - xe
        vec_ef(2) = ys - ye
        vec_ef(3) = zs - ze
        vec_ef(1:3) = vec_ef(1:3) / norm2(vec_ef(1:3))  ! normalize it

        ! Edge is convex
        if(dot_product(vec_ef(1:3), norm_1(1:3)) < 0.0 .and.  &
           dot_product(vec_ef(1:3), norm_2(1:3)) < 0.0) then
          edge_data(e) = 1

        ! Edge is concave
        else if(dot_product(vec_ef(1:3), norm_1(1:3)) > 0.0 .and.  &
                dot_product(vec_ef(1:3), norm_2(1:3)) > 0.0) then
          edge_data(e) = -1

        else
          print '(A,99F12.3)', 'BAD', norm2(vec_ef(1:3)), norm_1, norm_2,  &
                                      dot_product(norm_1(1:3), norm_2(1:3))
        end if

      end if

    end if  ! edge is between two boundary faces
    end if  ! edge is in this boundary condition

  end do  ! through edges


  !------------------------!
  !                        !
  !   Edge is on this BC   !
  !                        !
  !------------------------!
  do e = 1, Grid % n_edges
    if(Grid % edges_bc(bc, e) .gt. 0) then
    if( sum(Grid % edges_bc(1:Grid % n_bnd_regions, e)) .gt. 1 ) then
      if(edge_data(e) .eq. 0) then  ! hasn't been marked yet
        cnt = cnt + 1
        edge_data(e) = 1
      end if
    end if
    end if  ! it is in this bc
  end do

  ! Return a sum of all marked nodes
  N_Edges_In_Region = cnt

  end function

