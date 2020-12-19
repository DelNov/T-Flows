!==============================================================================!
  subroutine Grid_Mod_Calculate_Face_Geometry(grid)
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
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, pnt_to, pnt_from
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
  real    :: d_s, min_d, max_d
!==============================================================================!

  !----------------------------------------------!
  !   Calculate total surface of the cell face   !
  !----------------------------------------------!
  do s = 1, grid % n_faces
    grid % s(s) = sqrt(  grid % sx(s)*grid % sx(s)  &
                       + grid % sy(s)*grid % sy(s)  &
                       + grid % sz(s)*grid % sz(s) )
  end do

  !-------------------------------------------------------!
  !   Calculate the distance between neighbouring cells   !
  !    (For normal faces, including the periodic ones,    !
  !    dx, dy and dz are distanes between cell centers)   !
  !-------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    xc1 = grid % xc(c1)
    yc1 = grid % yc(c1)
    zc1 = grid % zc(c1)

    xc2 = grid % xc(c2) + grid % dx(s)
    yc2 = grid % yc(c2) + grid % dy(s)
    zc2 = grid % zc(c2) + grid % dz(s)

    grid % dx(s) = xc2-xc1
    grid % dy(s) = yc2-yc1
    grid % dz(s) = zc2-zc1
  end do  ! faces

  !---------------------------------------!
  !   Check distances stored in shadows   !
  !---------------------------------------!
  min_d = +HUGE
  max_d = -HUGE
  do s = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    d_s = sqrt(  grid % dx(s)*grid % dx(s)     &
               + grid % dy(s)*grid % dy(s)     &
               + grid % dz(s)*grid % dz(s) )
    min_d = min(min_d, d_s)
    max_d = max(max_d, d_s)
  end do
  call Comm_Mod_Global_Min_Real(min_d)
  call Comm_Mod_Global_Max_Real(max_d)
  if(this_proc < 2 .and. grid % n_shadows > 0) then
    print '(a,f9.3)', '# Minimum distance stored in shadow faces: ', min_d
    print '(a,f9.3)', '# Maximum distance stored in shadow faces: ', max_d
  end if

  !-----------------------------------------------------!
  !   Check distances between faces and their shadows   !
  !-----------------------------------------------------!
  min_d = +HUGE
  max_d = -HUGE
  do s = 1, grid % n_faces
    if(grid % faces_s(s) .ne. 0) then
      d_s = sqrt(  (grid % xf(s) - grid % xf(grid % faces_s(s)))**2     &
                 + (grid % yf(s) - grid % yf(grid % faces_s(s)))**2     &
                 + (grid % zf(s) - grid % zf(grid % faces_s(s)))**2 )
      min_d = min(min_d, d_s)
      max_d = max(max_d, d_s)
    end if
  end do
  call Comm_Mod_Global_Min_Real(min_d)
  call Comm_Mod_Global_Max_Real(max_d)
  if(this_proc < 2 .and. grid % n_shadows > 0) then
    print '(a,f9.3)', '# Minimum corrected distance at shadows:   ', min_d
    print '(a,f9.3)', '# Maximum corrected distance at shadows:   ', max_d
  end if

  !---------------------------------------------------!
  !   Counter and check pointer to and from shadows   !
  !---------------------------------------------------!
  pnt_to   = 0
  pnt_from = 0
  do s = 1, grid % n_faces
    if(grid % faces_s(s) .ne. 0) pnt_to = pnt_to + 1
  end do
  do s = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    if(grid % faces_s(s) .ne. 0) pnt_from = pnt_from + 1
  end do
  if(pnt_to .ne. pnt_from) then
    print *, '# Pointers to and from shadows wrong in processor: ', this_proc
    stop
  end if

  !--------------------------------------------!
  !   For shadows, dx, dy and dz are lengths   !
  !    of the periodic spans of the domain     !
  !--------------------------------------------!
  do s = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    xc1 = grid % xc(c1)
    yc1 = grid % yc(c1)
    zc1 = grid % zc(c1)

    xc2 = grid % xc(c2) + grid % dx(s)
    yc2 = grid % yc(c2) + grid % dy(s)
    zc2 = grid % zc(c2) + grid % dz(s)

    grid % dx(s) = grid % xf(grid % faces_s(s)) - grid % xf(s)
    grid % dy(s) = grid % yf(grid % faces_s(s)) - grid % yf(s)
    grid % dz(s) = grid % zf(grid % faces_s(s)) - grid % zf(s)
  end do  ! shadows

  !--------------------------------------------!
  !   Calculate total distance between cells   !
  !--------------------------------------------!
  do s = 1, grid % n_faces + grid % n_shadows
    grid % d(s) = sqrt(  grid % dx(s)*grid % dx(s)     &
                       + grid % dy(s)*grid % dy(s)     &
                       + grid % dz(s)*grid % dz(s) )
  end do

  !--------------------------------------------!
  !   Calculate weight factors for the faces   !
  !--------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Inside the flow, it has usual value: phi_f = f * phi_1 + (1-f) * phi_2
    grid % fw(s) = grid % f(s)

    ! Close to the wall, however, there is inversion. It takes
    ! the value from inside as the representative for the face.
    if(c2 < 0) then
      grid % fw(s) = 1.0
    end if
  end do

  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  grid % cell_near_wall(:) = .false.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        grid % cell_near_wall(c1) = .true.
      end if
    end if

  end do  ! faces

  end subroutine
