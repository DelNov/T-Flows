!==============================================================================!
  integer function N_Edges_On_Bnd_Color(grid, bc)
!------------------------------------------------------------------------------!
!   Counts and marks (with new_e) edges in the given boundary color            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bc
!-----------------------------------[Locals]-----------------------------------!
  integer :: e, cnt, s1, s2
  real    :: n1(3), n2(3)
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  grid % new_e(:) = 0

  do e = 1, grid % n_edges

    ! Edge is on this BC
    if(grid % edges_bc(bc, e) .gt. 0) then

      if(grid % edges_fb(1, e) .ne. 0 .and.  &
         grid % edges_fb(2, e) .ne. 0) then

        s1 = grid % edges_fb(1, e)
        s2 = grid % edges_fb(2, e)

        ! Compute normals of the two faces surrounding the edge
        n1(1) = grid % sx(s1);  n1(2) = grid % sy(s1);  n1(3) = grid % sz(s1)
        n2(1) = grid % sx(s2);  n2(2) = grid % sy(s2);  n2(3) = grid % sz(s2)
        n1(1:3) = n1(1:3) / norm2(n1(1:3))
        n2(1:3) = n2(1:3) / norm2(n2(1:3))

        ! PRINT '(6ES12.3)', DOT_PRODUCT(N1(1:3), N2(1:3))
        if(abs(dot_product(n1(1:3), n2(1:3))) < 0.7071) then
          cnt = cnt + 1
          grid % new_e(e) = cnt
        end if
      end if

    end if

    ! Edge is on this BC
    if(grid % edges_bc(bc, e) .gt. 0) then

      ! But there is more
      if( sum(grid % edges_bc(1:grid % n_bnd_cond, e)) .gt. 1 ) then
        cnt = cnt + 1
        grid % new_e(e) = cnt
      end if
    end if
  end do

  ! Return a sum of all marked nodes
  N_Edges_On_Bnd_Color = cnt

  end function

