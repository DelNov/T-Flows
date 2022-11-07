!==============================================================================!
  subroutine Mark_Cells_And_Faces(Front, phi)
!------------------------------------------------------------------------------!
!   Mark cells and find faces at surface                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Var_Type),    target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  integer                   :: c1, c2, e, s, n_fac
  real                      :: phi_c1, phi_c2, w1, w2
! character(SL)             :: line1, line2, line3
!==============================================================================!

  ! Take aliases
  Flow  => Front % pnt_flow
  Grid  => Front % pnt_grid

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  Front % elem_in_cell(:) = 0  ! not at surface
  do e = 1, Front % n_elems
    Front % elem_in_cell(Front % Elem(e) % cell) = e
  end do

  !-----------!
  !           !
  !   Faces   !
  !           !
  !-----------!
  n_fac = 0
  Front % elems_at_face(:,:) = 0  ! not at surface
  Grid % xs(:) = 0.0
  Grid % ys(:) = 0.0
  Grid % zs(:) = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_c1 = phi % n(c1)
    phi_c2 = phi % n(c2)

    if(Math % Approx_Real(phi_c1, 0.5)) phi_c1 = phi_c1 - MICRO
    if(Math % Approx_Real(phi_c2, 0.5)) phi_c2 = phi_c2 - MICRO

    !---------------------------------------!
    !   If face crosses the "phi_e" value   !
    !---------------------------------------!
    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then

      ! There is a ludicrously simple way to find the coordintate of the
      ! intersection between surrounding cell connections and element:
      ! By a linear interpolation :-/
      w1 = abs(phi_c2 - 0.5) / abs(phi_c1 - phi_c2)
      w2 = abs(phi_c1 - 0.5) / abs(phi_c1 - phi_c2)

      ! Intersection point (dx, dy and dz are to take care of periodicity)
      Grid % xs(s) = w1 * Grid % xc(c1) + w2 * (Grid % xc(c1) + Grid % dx(s))
      Grid % ys(s) = w1 * Grid % yc(c1) + w2 * (Grid % yc(c1) + Grid % dy(s))
      Grid % zs(s) = w1 * Grid % zc(c1) + w2 * (Grid % zc(c1) + Grid % dz(s))

      ! Store the elements at this face: essentially just
      ! copy the elements at cells surrounding this face
      n_fac = n_fac + 1
      Front % elems_at_face(1, s) = Front % elem_in_cell(c1)
      Front % elems_at_face(2, s) = Front % elem_in_cell(c2)

    end if  ! face crosses 0.5
  end do    ! through faces

  ! This check probably doesn't make a hell of a lot of sense
  ! if(Front % n_elems .ne. n_fac) then
  !   write(line1, '(a)')    " It seems that not all face "  //  &
  !                          " intersections have been found! \n "
  !   write(line2, '(a,i6)') " Front % n_elems = ", Front % n_elems
  !   write(line3, '(a,i6)') " n_fac           = ", n_fac
  !   call Message % Warning(60, line1 // line2 // line3,  &
  !                          file=__FILE__, line=__LINE__)
  !
  ! end if

  end subroutine
