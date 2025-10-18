!==============================================================================!
  subroutine Calculate_Face_Geometry(Grid)
!------------------------------------------------------------------------------!
!>  Calculates additional face-base geometrical quantities for Process.
!------------------------------------------------------------------------------!
!   This subroutine is called only from Process, and it was originally used    !
!   whenever I was lazy to change the format of .cfn (former .cns) files.      !
!   However, in the meanwhile I realized it is not so bad after all.  There    !
!   are quantities which can be defined only in processors, for example when   !
!   something depends on boundary conditions.                                  !
!                                                                              !
!   More importantly, certaing integrity checks can be done only here!         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, i_cel, i_fac, sh, n, pnt_to, pnt_from, fail_count
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
  real    :: d_s, min_d, max_d, sx, sy, sz
  real    :: dist(3), surf(3)
!==============================================================================!

  if(First_Proc()) print '(a)', ' # Checking the integrity of cell faces ...'

  !-----------------------------------!
  !   Perform some integrity checks   !
  !-----------------------------------!

  ! Check integrity of internal face surfaces
  fail_count = 0
  do s = 1, Grid % n_faces
    call Grid % Faces_Surface(s, sx, sy, sz)
    if(Grid % sx(s) * sx + Grid % sy(s) * sy + Grid % sz(s) * sz < 0.0) then
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      fail_count = fail_count + 1
    end if
  end do
  if(fail_count .gt. 0) call Message % Error(60,                               &
                            'Face integrity test for internal faces failed!',  &
                             file=__FILE__)

  ! Check integrity of shadow face surfaces
  fail_count = 0
  do sh = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    call Grid % Faces_Surface(sh, sx, sy, sz)
    if(Grid % sx(sh) * sx + Grid % sy(sh) * sy + Grid % sz(sh) * sz < 0.0) then
      fail_count = fail_count + 1
    end if
  end do
  if(fail_count .gt. 0) call Message % Error(60,                             &
                            'Face integrity test for shadow faces failed!',  &
                             file=__FILE__)

  ! Check integrity of shadow faces connectivity (1)
  fail_count = 0
  do s = 1, Grid % n_faces
    sh = Grid % faces_s(s)
    if(sh .ne. 0) then
      if(sh .le. Grid % n_faces) fail_count = fail_count + 1
      if(Grid % faces_c(1,s) .ne. Grid % faces_c(1,sh) .or.  &
         Grid % faces_c(2,s) .ne. Grid % faces_c(2,sh))      &
        fail_count = fail_count + 1
    end if
  end do
  if(fail_count .gt. 0) call Message % Error(60,                        &
                            'Shadow face connectivity test 1 failed!',  &
                             file=__FILE__)

  ! Check integrity of shadow faces connectivity (2)
  fail_count = 0
  do s = 1, Grid % n_faces
    sh = Grid % faces_s(s)
    if(sh .ne. 0) then
      if(Grid % faces_n_nodes(s) .ne. Grid % faces_n_nodes(sh))  &
        fail_count = fail_count + 1
      n = Grid % faces_n_nodes(sh)
      if(any(Grid % faces_n(1:n, s) .eq. Grid % faces_n(1:n,sh)))  &
        fail_count = fail_count + 1
      if(any(Grid % faces_n(1:n, sh) .le. 0))  &
        fail_count = fail_count + 1
    end if
  end do
  if(fail_count .gt. 0) call Message % Error(60,                        &
                            'Shadow face connectivity test 2 failed!',  &
                             file=__FILE__)

  !------------------------------------------!
  !                                          !
  !   Check integrity of cells' enclosures   !
  !                                          !
  !------------------------------------------!
  if(DEBUG) call Grid % Check_Cells_Closure()

  !----------------------------------------------!
  !                                              !
  !   Calculate total surface of the cell face   !
  !                                              !
  !----------------------------------------------!
  do s = 1, Grid % n_faces
    Grid % s(s) = sqrt(  Grid % sx(s)*Grid % sx(s)  &
                       + Grid % sy(s)*Grid % sy(s)  &
                       + Grid % sz(s)*Grid % sz(s) )
  end do

  !-------------------------------------------------------!
  !                                                       !
  !   Calculate the distance between neighbouring cells   !
  !    (For normal faces, including the periodic ones,    !
  !    dx, dy and dz are distanes between cell centers)   !
  !                                                       !
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

  ! Check integrity of faces surface and distance vectors
  fail_count = 0
  do s = 1, Grid % n_faces
    dist(1) = Grid % dx(s)
    dist(2) = Grid % dy(s)
    dist(3) = Grid % dz(s)
    surf(1) = Grid % sx(s)
    surf(2) = Grid % sy(s)
    surf(3) = Grid % sz(s)

    ! Check the dot product of face surface and cell connection along the way
    if(dot_product(dist, surf) < 0.0) then
      fail_count = fail_count + 1
    end if

  end do

  if(fail_count .gt. 0) call Message % Error(40,                            &
                            "The face or cell connection orientation "  //  &
                            "is wrong. The code can't continue like  "  //  &
                            "this and will terminate now.",                 &
                             file=__FILE__, line=__LINE__)

  if(First_Proc()) print '(a)', ' # All integrity tests passed'

  !---------------------------------------!
  !                                       !
  !   Check distances stored in shadows   !
  !                                       !
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
  call Global % Min_Real(min_d)
  call Global % Max_Real(max_d)
  if(First_Proc() .and. Grid % n_shadows > 0) then
    print '(a,f9.3)', ' # Minimum distance stored in shadow faces: ', min_d
    print '(a,f9.3)', ' # Maximum distance stored in shadow faces: ', max_d
  end if

  !---------------------------------------------------------!
  !                                                         !
  !   Set up straight boundary conditions for periodicity   !
  !   (This is important for copy boundary conditions, do   !
  !        not erase this thinking it is not needed)        !
  !                                                         !
  !---------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(Grid % faces_s(s) .ne. 0) then
      sh = Grid % faces_s(s)

      if(abs(Grid % dx(s)) > PICO) then
        if(abs(Grid % xc(c2) - Grid % xc(c1) ) > 1.5 * abs(Grid % dx(s))) then
          Assert(Grid % faces_s(s) > 0)
          Grid % region % at_face(s) = Grid % per_x_reg
        end if
      end if

      if(abs(Grid % dy(s)) > PICO) then
        if(abs(Grid % yc(c2) - Grid % yc(c1) ) > 1.5 * abs(Grid % dy(s))) then
          Assert(Grid % faces_s(s) > 0)
          Grid % region % at_face(s) = Grid % per_y_reg
        end if
      end if

      if(abs(Grid % dz(s)) > PICO) then
        if(abs(Grid % zc(c2) - Grid % zc(c1) ) > 1.5 * abs(Grid % dz(s))) then
          Assert(Grid % faces_s(s) > 0)
          Grid % region % at_face(s) = Grid % per_z_reg
        end if
      end if

    end if
  end do

  !-----------------------------------------------------!
  !                                                     !
  !   Check distances between faces and their shadows   !
  !                                                     !
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
  call Global % Min_Real(min_d)
  call Global % Max_Real(max_d)
  if(First_Proc() .and. Grid % n_shadows > 0) then
    print '(a,f9.3)', ' # Minimum corrected distance at shadows:   ', min_d
    print '(a,f9.3)', ' # Maximum corrected distance at shadows:   ', max_d
  end if

  !---------------------------------------------------!
  !                                                   !
  !   Counter and check pointer to and from shadows   !
  !                                                   !
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
    print *, '# Pointers to and from shadows wrong in processor: ', This_Proc()
    stop
  end if

  !--------------------------------------------!
  !                                            !
  !   For shadows, dx, dy and dz are lengths   !
  !    of the periodic spans of the domain     !
  !                                            !
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
  !                                            !
  !   Calculate total distance between cells   !
  !                                            !
  !--------------------------------------------!
  do s = 1, Grid % n_faces + Grid % n_shadows
    Grid % d(s) = sqrt(  Grid % dx(s)*Grid % dx(s)     &
                       + Grid % dy(s)*Grid % dy(s)     &
                       + Grid % dz(s)*Grid % dz(s) )
  end do

  !--------------------------------------------!
  !                                            !
  !   Calculate weight factors for the faces   !
  !                                            !
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
  !                              !
  !   Find the near-wall cells   !
  !                              !
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

  !--------------------------------------------------!
  !                                                  !
  !   Fill up the cells_n_cells and cells_c arrays   !
  !                                                  !
  !--------------------------------------------------!
  Grid % cells_n_cells(:) = 0
  do c1 = -Grid % n_bnd_cells, Grid % n_cells

    if(c1 .eq. 0) cycle  ! skip cell zero

    n = 0  ! counter for Grid % cells_n_cells
    do i_fac = 1, Grid % cells_n_faces(c1)
      s = Grid % cells_f(i_fac, c1)
      n = n + 1

      call Enlarge % Matrix_Int(Grid % cells_c, i = (/1, n/))
      c2 = Grid % faces_c(1,s) + Grid % faces_c(2,s) - c1
      Grid % cells_c(n, c1) = c2
    end do
    Grid % cells_n_cells(c1) = n

    ! Sort cell's neighbours and carry faces along
    call Sort % Int_Carry_Int(Grid % cells_c(1:n, c1),  &
                              Grid % cells_f(1:n, c1))

    ! Find the first "inside" cell
    Grid % cells_i_cells(c1) = 0
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      if(c2 .gt. 0) then
        Grid % cells_i_cells(c1) = i_cel
        exit
      end if
    end do
    Assert(Grid % cells_i_cells(c1) .ne. 0)

  end do

  !-------------------------------------!
  !   Check the Grid % cells_c matrix   !
  !-------------------------------------!

  ! Test 1
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % Comm % cell_proc(c) .eq. This_Proc()) then
      Assert(Grid % cells_n_cells(c) .eq. Grid % cells_n_faces(c))
    end if
  end do

  ! Test 2
  do c1 = -Grid % n_bnd_cells, Grid % n_cells
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      Assert(Grid % faces_c(1,s) + Grid % faces_c(2,s) .eq. c1 + c2)
    end do
  end do

  ! Test 3
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do i_cel = 1, Grid % cells_n_cells(c)
      Assert(Grid % cells_c(i_cel, c) .ne. 0)
    end do
  end do

  end subroutine
