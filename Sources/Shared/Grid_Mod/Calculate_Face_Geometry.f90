!==============================================================================!
  subroutine Grid_Mod_Calculate_Face_Geometry(grid)
!------------------------------------------------------------------------------!
!   Calculates additional face-base geometrical quantities for Process.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
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
