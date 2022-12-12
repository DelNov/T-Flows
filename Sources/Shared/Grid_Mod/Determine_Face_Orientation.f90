!==============================================================================!
  subroutine Determine_Face_Orientation(Grid)
!------------------------------------------------------------------------------!
!   Face orientation is important only for Isoap, hence for VOF simulations    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, n1, n2, i_fac
  real    :: vec1(3), vec2(3), prod(3), dist(3), surf(3)
!==============================================================================!

  !--------------------------------------!
  !   Check orientaion of faces' nodes   !
  !--------------------------------------+------------------------------!
  !   It is quite clear that faces's nodes could be oriented in both    !
  !   clockwise and counter-clockwise direction, depending on the       !
  !   software which was used to create the mesh, all the algorithms    !
  !   which have been used to extract faces, create polyhedral grids,   !
  !   so maybe it is the best to check them all here, in Process.       !
  !---------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    dist(1) = Grid % dx(s)
    dist(2) = Grid % dy(s)
    dist(3) = Grid % dz(s)
    surf(1) = Grid % sx(s)
    surf(2) = Grid % sy(s)
    surf(3) = Grid % sz(s)

    ! Check the dot product of face surface and cell connection along the way
    if(dot_product(dist, surf) < 0.0) then
      call Message % Error(40,                                  &
                "The face or cell connection orientation "  //  &
                "is wrong. The code can't continue like  "  //  &
                "this and will terminate now.",                 &
                file=__FILE__, line=__LINE__)
    end if

    ! Browsing through all the nodes would clearly be
    ! an overkill here, so we take only the first two
    n1 = Grid % faces_n(1, s)
    n2 = Grid % faces_n(2, s)
    vec1(1) = Grid % xn(n1) - Grid % xf(s)
    vec1(2) = Grid % yn(n1) - Grid % yf(s)
    vec1(3) = Grid % zn(n1) - Grid % zf(s)
    vec2(1) = Grid % xn(n2) - Grid % xf(s)
    vec2(2) = Grid % yn(n2) - Grid % yf(s)
    vec2(3) = Grid % zn(n2) - Grid % zf(s)
    prod = Math % Cross_Product(vec1, vec2)
    if(dot_product(prod, dist) < 0.0) then
      call Sort % Reverse_Order_Int(Grid % faces_n(1:Grid % faces_n_nodes(s),s))
    end if
  end do

  !-----------------------------------------------------------------!
  !   Find if faces around cells are oriented inwards or outwards   !
  !-----------------------------------------------------------------!
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells

    ! Even Grid % cells_f_orient are initially allocated with very low dimension
    call Adjust_First_Dim(Grid % cells_n_faces(c), Grid % cells_f_orient)

    do i_fac = 1, Grid % cells_n_faces(c)  ! browse through faces of the cell
      s  = Grid % cells_f(i_fac, c)        ! take face's true number
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! Skip the shadow faces, no reliable
      ! information on c1 and c2 on them
      if(s .lt. Grid % n_faces - Grid % n_shadows) then

        ! Face orientation (see above) and surface vector point
        ! from c1 to c2.  So, if c1 is equal to c, it means that
        ! the relative surface i_fac to c is pointing outwards
        if(c1 .eq. c) then
          Grid % cells_f_orient(i_fac, c) = OUTWARDS
        else if(c2 .eq. c) then
          Grid % cells_f_orient(i_fac, c) = INWARDS
        else
          call Message % Error(40,                          &
                  "Something is seriously wrong with "  //  &
                  "face to cell connectivity.",             &
                  file=__FILE__, line=__LINE__)
        end if

      end if  ! face is not a shadow face

    end do  ! through cell's face
  end do    ! through cells

  end subroutine
