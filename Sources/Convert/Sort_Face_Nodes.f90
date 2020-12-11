!==============================================================================!
  subroutine Sort_Face_Nodes(grid, s, concave_link)
!------------------------------------------------------------------------------!
!   Sort nodes of a given face in rotational fashion                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: s
  integer         :: concave_link(2, grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n, nn, cnt
  real    :: normal_p(3), center_p(3), x_p(3), y_p(3), sense(3)
  real    :: rp_3d(3,MAX_FACES_N_NODES), np_3d(3,MAX_FACES_N_NODES)
  real    :: rp_2d(2,MAX_FACES_N_NODES)
  real    :: sorting(MAX_FACES_N_NODES)  ! sorting criterion
  integer :: order(MAX_FACES_N_NODES)    ! carry-on array with indices
  integer :: max_loc(2), conc_n
  real    :: conc_xn, conc_yn, conc_zn
  real    :: prod(3), prod_mag, angles(MAX_FACES_N_NODES, MAX_FACES_N_NODES)
  real    :: sumang, criter
! integer :: k, ni, nj, nk, min_loc
! real    :: vec_ji(3), vec_jk(3), mag_ji, mag_jk, dot_prod
!==============================================================================!

  ! Take alias
  nn = grid % faces_n_nodes(s)

  !---------------------------------------------------------------!
  !   A way to deal with concave faces.  Mark its concave node,   !
  !   move it to a position which makes the face convex, use      !
  !   the sorting algorithm developed before, and then move the   !
  !   concave node back to its original position.  This is good   !
  !   only for faces wiht one concave node.  If it has two, it    !
  !   is likely that the entire face is in an edge, between two   !
  !   planes.  This is step 1.  Step 2 is after sorting.          !
  !---------------------------------------------------------------!
  cnt    = 0  ! not really needed; it seems that each face can have only one
  conc_n = 0
  do i = 1, nn
    n = grid % faces_n(i, s)
    if(concave_link(1, n) .gt. 0 .and. concave_link(2, n) .gt. 0) then
      conc_n = n
      cnt    = cnt + 1

      ! Store concave node coordinates
      conc_xn = grid % xn(n)
      conc_yn = grid % yn(n)
      conc_zn = grid % zn(n)

      ! Move the concave node in a convex position
      grid % xn(n) = 0.5 * (  grid % xn(concave_link(1, n))     &
                            + grid % xn(concave_link(2, n))  )
      grid % yn(n) = 0.5 * (  grid % yn(concave_link(1, n))     &
                            + grid % yn(concave_link(2, n))  )
      grid % zn(n) = 0.5 * (  grid % zn(concave_link(1, n))     &
                            + grid % zn(concave_link(2, n))  )
    end if
  end do
  ! if(cnt .eq. 1) print *, '# Found a concave link in face', s
  ! if(cnt .eq. 2) print *, '# Found two concave links in face', s

  !-----------------------------------!
  !   Find the plane's center point   !
  !-----------------------------------!
  center_p(:) = 0
  do i = 1, nn
    n = grid % faces_n(i, s)
    center_p(1) = center_p(1) + grid % xn(n)
    center_p(2) = center_p(2) + grid % yn(n)
    center_p(3) = center_p(3) + grid % zn(n)
  end do
  center_p(1:3) = center_p(1:3) / real(nn)

  !------------------------------------------------------------------!
  !   Find nodes' relative positions to the center just calculated   !
  !------------------------------------------------------------------!
  do i = 1, nn  ! use "i", not "i_nod" here
    n = grid % faces_n(i, s)
    rp_3d(1, i) = grid % xn(n) - center_p(1)
    rp_3d(2, i) = grid % yn(n) - center_p(2)
    rp_3d(3, i) = grid % zn(n) - center_p(3)
  end do

  !--------------------------------!
  !   Find a normal of the plane   !
  !--------------------------------!

  ! Calculate normalized relative vectors
  do i = 1, nn  ! use "i", not "i_nod" here
    np_3d(1:3,i) = rp_3d(1:3,i) / norm2(rp_3d(1:3,i))
  end do

  ! Multiply each with each to find angles between them
  ! (Each of these multiplications is a guess for the normal)
  cnt               = 0
  sumang            = 0.0
  angles(1:nn,1:nn) = 0.0
  do i = 1, nn
    do j = i + 1, nn
      prod(1:3) = Math_Mod_Cross_Product(np_3d(1:3,i), np_3d(1:3,j))
      prod_mag = min(norm2(prod(1:3)), 1.0-MILI)
      angles(i,j) = asin(prod_mag) * 57.2957795131
      if(angles(i,j) > MILI) then
        cnt = cnt + 1
        sumang = sumang + angles(i,j)
      end if
    end do
  end do
  criter = 0.75 * sumang / real(cnt)

  ! Find maximum angle, that one wil be relevant for the sense od normal
  max_loc = maxloc(angles(1:nn,1:nn));  i = max_loc(1);  j = max_loc(2)
  sense(:) = Math_Mod_Cross_Product(np_3d(1:3,i), np_3d(1:3,j))

  ! Multiply each with each selectivelly, taking into account only bigger
  ! angles, and taking care to correct the signs in proper sense.
  ! Average the surface normal along the way.
  cnt           = 0
  normal_p(1:3) = 0.0
  do i = 1, nn
    do j = i + 1, nn
      if(angles(i,j) > criter) then  ! take only reasonalby big angles
        cnt = cnt + 1                ! one more sample
        prod(1:3) = Math_Mod_Cross_Product(np_3d(1:3,i), np_3d(1:3,j))
        if(dot_product(prod(1:3), sense(1:3)) < 0) then  ! correct the sign ...
          prod(1:3) = -prod(1:3)                         ! ... if needed
        end if
        normal_p(1:3) = normal_p(1:3) + prod(1:3)
      end if
    end do
  end do
  normal_p(1:3) = normal_p(1:3) / real(cnt)

  !------------------------------------------------!
  !   Create a 2D coordinate system on the plane   !
  !   using its central point and normal, and      !
  !   then project all the points on that system   !
  !------------------------------------------------!

  ! Define x-axis in the plane
  x_p(1:3) = Math_Mod_Cross_Product(normal_p(1:3), rp_3d(1:3, 1))
  x_p(1:3) = -x_p(1:3) / norm2(x_p(1:3))

  ! Define y-axis in the plane
  y_p(1:3) = Math_Mod_Cross_Product(normal_p(1:3), x_p(1:3))
  y_p(1:3) = y_p(1:3) / norm2(y_p(1:3))

  ! Project relative points to the plane's coordinate system
  do i = 1, nn  ! use "i", not "i_nod" here
    n = grid % faces_n(i, s)
    rp_2d(1, i) = dot_product(x_p(1:3), rp_3d(1:3, i))
    rp_2d(2, i) = dot_product(y_p(1:3), rp_3d(1:3, i))
  end do

  do i = 1, nn  ! use "i", not "i_nod" here
    sorting(i) = atan2(rp_2d(1,i), rp_2d(2,i)) * 57.2957795131
    order(i) = grid % faces_n(i, s)
  end do

  call Sort_Mod_Real_Carry_2_Int(sorting(1:nn),           &
                                 grid % faces_n(1:nn,s),  &
                                 order(1:nn))

  !---------------------------------------------------------------!
  !   A way to deal with concave faces.  Mark its concave node,   !
  !   move it to a position which makes the face convex, use      !
  !   the sorting algorithm developed before, and then move the   !
  !   concave node back to its original position.  This is good   !
  !   only for faces wiht one concave node.  If it has two, it    !
  !   is likely that the entire face is in an edge, between two   !
  !   planes.  This is step 2.  Step 1 is at the top.             !
  !---------------------------------------------------------------!
  if(conc_n .ne. 0) then
    do i = 1, nn
      n = grid % faces_n(i, s)
      if(n .eq. conc_n) then
        grid % xn(n) = conc_xn
        grid % yn(n) = conc_yn
        grid % zn(n) = conc_zn
      end if
    end do
  end if

  ! Rotate for better visualisation with Paraview.  Paraview constructs
  ! a polygon from triangles which all originate in the first node.  If
  ! the first node is defined exactly at the sharpest edge, which is also
  ! the first one to be inserted in polynomial face in the calling func.
  ! All in all: try to keep the first node first for better visuals.
  if(conc_n .ne. 0) then
    do while(order(1) .ne. conc_n)
      order(1:nn)             = cshift(order(1:nn), 1)
      grid % faces_n(1:nn, s) = cshift(grid % faces_n(1:nn, s), 1)
    end do
  end if

  end subroutine
