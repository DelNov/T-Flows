!==============================================================================!
  integer function N_Sharp_Edges(grid)
!------------------------------------------------------------------------------!
!   Counts and marks (with new_e) edges at sharp boundaries and in between     !
!   different boundary conditions.                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: e, cnt, s1, s2
  real    :: norm_1(3), norm_2(3)
!==============================================================================!

  ! Nullify on entry
  cnt             = 0
  grid % new_e(:) = 0

  !---------------------------------------------!
  !   Fetch geometrically sharp edges (first)   !
  !---------------------------------------------!
  do e = 1, grid % n_edges

    ! If edge is in-between two boundary faces
    ! (yes, grid % edges_f stores only boundary faces, very bad name choice)
    if(grid % edges_fb(1, e) .ne. 0 .and.  &
       grid % edges_fb(2, e) .ne. 0) then

      s1 = grid % edges_fb(1, e)
      s2 = grid % edges_fb(2, e)

      ! Compute normals of the two faces surrounding the edge
      norm_1(1) = grid % sx(s1);  norm_2(1) = grid % sx(s2)
      norm_1(2) = grid % sy(s1);  norm_2(2) = grid % sy(s2)
      norm_1(3) = grid % sz(s1);  norm_2(3) = grid % sz(s2)
      norm_1(1:3) = norm_1(1:3) / norm2(norm_1(1:3))
      norm_2(1:3) = norm_2(1:3) / norm2(norm_2(1:3))

      ! This has to do with 45 degrees, but I forgot how exactly
      if(abs(dot_product(norm_1(1:3), norm_2(1:3))) < 0.7071) then
        cnt = cnt + 1
        grid % new_e(e) = cnt
      end if
    end if

  end do

  !-------------------------------------------------------------------!
  !   But also fetch edges in-between different boundary conditions   !
  !-------------------------------------------------------------------!
  do e = 1, grid % n_edges
    if( sum(grid % edges_bc(1:grid % n_bnd_cond, e)) .gt. 1 ) then
      if(grid % new_e(e) .eq. 0) then  ! hasn't been marked yet
        cnt = cnt + 1
        grid % new_e(e) = cnt
      end if
    end if
  end do

  ! Return a sum of all marked nodes
  N_Sharp_Edges = cnt

  end function

