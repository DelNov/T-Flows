!==============================================================================!
  subroutine Calculate_Face_Geometry(Grid)
!------------------------------------------------------------------------------!
!   Calculates additional face-base geometrical quantities for Process.        !
!                                                                              !
!   This subroutine is called only from Process, and it was originally used    !
!   whenever I was lazy to change the format of .cfn (former .cns) files.      !
!   That is very unfortunate, and I should use it as little as possible and    !
!   never expand it.  Actually, I should try to get rid of it completelly.     !
!                                                                              !
!   I am not quite sure, but it is conceiveable that shadow faces were         !
!   introduced only to make sure this routine works as it should. Bad!         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, sh, pnt_to, pnt_from, i, n1, n2, i_fac
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
  real    :: d_s, min_d, max_d
  real    :: vec1(3), vec2(3), prod(3), dist(3), surf(3)
!==============================================================================!

  !----------------------------------------------!
  !   Calculate total surface of the cell face   !
  !----------------------------------------------!
  do s = 1, Grid % n_faces
    Grid % s(s) = sqrt(  Grid % sx(s)*Grid % sx(s)  &
                       + Grid % sy(s)*Grid % sy(s)  &
                       + Grid % sz(s)*Grid % sz(s) )
  end do

  !-------------------------------------------------------!
  !   Calculate the distance between neighbouring cells   !
  !    (For normal faces, including the periodic ones,    !
  !    dx, dy and dz are distanes between cell centers)   !
  !-------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    xc1 = Grid % xc(c1)
    yc1 = Grid % yc(c1)
    zc1 = Grid % zc(c1)

    xc2 = Grid % xc(c2) + Grid % dx(s)
    yc2 = Grid % yc(c2) + Grid % dy(s)
    zc2 = Grid % zc(c2) + Grid % dz(s)

    Grid % dx(s) = xc2-xc1
    Grid % dy(s) = yc2-yc1
    Grid % dz(s) = zc2-zc1
  end do  ! faces

  !---------------------------------------!
  !   Check distances stored in shadows   !
  !---------------------------------------!
  min_d = +HUGE
  max_d = -HUGE
  do s = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    d_s = sqrt(  Grid % dx(s)*Grid % dx(s)     &
               + Grid % dy(s)*Grid % dy(s)     &
               + Grid % dz(s)*Grid % dz(s) )
    min_d = min(min_d, d_s)
    max_d = max(max_d, d_s)
  end do
  call Comm_Mod_Global_Min_Real(min_d)
  call Comm_Mod_Global_Max_Real(max_d)
  if(this_proc < 2 .and. Grid % n_shadows > 0) then
    print '(a,f9.3)', ' # Minimum distance stored in shadow faces: ', min_d
    print '(a,f9.3)', ' # Maximum distance stored in shadow faces: ', max_d
  end if

  !---------------------------------------------------------!
  !   Set up straight boundary conditions for periodicity   !
  !   (This is important for copy boundary conditions, do   !
  !        not erase this thinking it is not needed)        !
  !---------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(Grid % faces_s(s) .ne. 0) then
      sh = Grid % faces_s(s)

      if(abs(Grid % xc(c2) - Grid % xc(c1) ) > 1.5 * abs(Grid % dx(s))) then
        if(Grid % xc(c2) > Grid % xc(c1))                      &
          Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 1
      end if

      if(abs(Grid % yc(c2) - Grid % yc(c1) ) > 1.5 * abs(Grid % dy(s))) then
        if(Grid % yc(c2) > Grid % yc(c1))                      &
          Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 2
      end if

      if(abs(Grid % zc(c2) - Grid % zc(c1) ) > 1.5 * abs(Grid % dz(s))) then
        if(Grid % zc(c2) > Grid % zc(c1))                      &
          Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 3
      end if

    end if
  end do

  !-----------------------------------------------------!
  !   Check distances between faces and their shadows   !
  !-----------------------------------------------------!
  min_d = +HUGE
  max_d = -HUGE
  do s = 1, Grid % n_faces
    if(Grid % faces_s(s) .ne. 0) then
      d_s = sqrt(  (Grid % xf(s) - Grid % xf(Grid % faces_s(s)))**2     &
                 + (Grid % yf(s) - Grid % yf(Grid % faces_s(s)))**2     &
                 + (Grid % zf(s) - Grid % zf(Grid % faces_s(s)))**2 )
      min_d = min(min_d, d_s)
      max_d = max(max_d, d_s)
    end if
  end do
  call Comm_Mod_Global_Min_Real(min_d)
  call Comm_Mod_Global_Max_Real(max_d)
  if(this_proc < 2 .and. Grid % n_shadows > 0) then
    print '(a,f9.3)', ' # Minimum corrected distance at shadows:   ', min_d
    print '(a,f9.3)', ' # Maximum corrected distance at shadows:   ', max_d
  end if

  !---------------------------------------------------!
  !   Counter and check pointer to and from shadows   !
  !---------------------------------------------------!
  pnt_to   = 0
  pnt_from = 0
  do s = 1, Grid % n_faces
    if(Grid % faces_s(s) .ne. 0) pnt_to = pnt_to + 1
  end do
  do s = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    if(Grid % faces_s(s) .ne. 0) pnt_from = pnt_from + 1
  end do
  if(pnt_to .ne. pnt_from) then
    print *, '# Pointers to and from shadows wrong in processor: ', this_proc
    stop
  end if

  !--------------------------------------------!
  !   For shadows, dx, dy and dz are lengths   !
  !    of the periodic spans of the domain     !
  !--------------------------------------------!
  do s = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    xc1 = Grid % xc(c1)
    yc1 = Grid % yc(c1)
    zc1 = Grid % zc(c1)

    xc2 = Grid % xc(c2) + Grid % dx(s)
    yc2 = Grid % yc(c2) + Grid % dy(s)
    zc2 = Grid % zc(c2) + Grid % dz(s)

    Grid % dx(s) = Grid % xf(Grid % faces_s(s)) - Grid % xf(s)
    Grid % dy(s) = Grid % yf(Grid % faces_s(s)) - Grid % yf(s)
    Grid % dz(s) = Grid % zf(Grid % faces_s(s)) - Grid % zf(s)
  end do  ! shadows

  !--------------------------------------------!
  !   Calculate total distance between cells   !
  !--------------------------------------------!
  do s = 1, Grid % n_faces + Grid % n_shadows
    Grid % d(s) = sqrt(  Grid % dx(s)*Grid % dx(s)     &
                       + Grid % dy(s)*Grid % dy(s)     &
                       + Grid % dz(s)*Grid % dz(s) )
  end do

  !--------------------------------------------!
  !   Calculate weight factors for the faces   !
  !--------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Inside the flow, it has usual value: phi_f = f * phi_1 + (1-f) * phi_2
    Grid % fw(s) = Grid % f(s)

    ! Close to the wall, however, there is inversion. It takes
    ! the value from inside as the representative for the face.
    if(c2 < 0) then
      Grid % fw(s) = 1.0
    end if
  end do

  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  Grid % cell_near_wall(:) = .false.

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if(Bnd_Cond_Type(Grid,c2) .eq. WALL .or.  &
         Bnd_Cond_Type(Grid,c2) .eq. WALLFL) then
        Grid % cell_near_wall(c1) = .true.
      end if
    end if

  end do  ! faces

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
