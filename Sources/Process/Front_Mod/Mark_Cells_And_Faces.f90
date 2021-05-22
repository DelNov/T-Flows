!==============================================================================!
  subroutine Mark_Cells_And_Faces(Front, phi)
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Var_Type),    target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  integer                   :: c, c1, c2, i_cel, e, s
  integer                   :: i_ver, j_ver, i, j, n_fac, n_int
  real                      :: nx, ny, nz, l, lx, ly, lz, dsc1, phi_c1, phi_c2
  real                      :: vec_i(3), vec_j(3), vec_ixj(3), xs, ys, zs
!==============================================================================!

  ! Take aliases
  flow  => Front % pnt_flow
  grid  => Front % pnt_grid

  !------------------------------------------!
  !                                          !
  !   Mark cells and find faces at surface   !
  !                                          !
  !------------------------------------------!
  if(flow % mass_transfer) then

    !-----------!
    !   Cells   !
    !-----------!
    Front % cell_at_elem(:) = 0  ! not at surface
    do e = 1, Front % n_elems
      Front % cell_at_elem(Front % elem(e) % cell) = e
    end do

    !-----------!
    !   Faces   !
    !-----------!
    n_fac = 0
    Front % face_at_elem(:,:) = 0  ! not at surface
    grid % xs(:) = 0.0
    grid % ys(:) = 0.0
    grid % zs(:) = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      phi_c1 = phi % n(c1)
      phi_c2 = phi % n(c2)

      if(Math_Mod_Approx_Real(phi_c1, 0.5)) phi_c1 = phi_c1 - MICRO
      if(Math_Mod_Approx_Real(phi_c2, 0.5)) phi_c2 = phi_c2 - MICRO

      !---------------------------------------!
      !   If face crosses the "phi_e" value   !
      !---------------------------------------!
      if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then

        ! Vector connecting cell centers c1 and c2
        lx = grid % dx(s)
        ly = grid % dy(s)
        lz = grid % dz(s)
        l  = sqrt(lx**2 + ly**2 + lz**2)

        n_int = 0
        do i_cel = 1, 2
          c = grid % faces_c(i_cel, s)

          ! Cell c contains surface
          e = Front % cell_at_elem(c)
          if(e .ne. 0) then
            nx = Front % elem(e) % nx
            ny = Front % elem(e) % ny
            nz = Front % elem(e) % nz

            ! If cell contains surface and angle is not bigger than asin(15)
            if( (lx*nx + ly*ny + lz*nz)/l > 0.258819 ) then

              nx = Front % elem(e) % nx
              ny = Front % elem(e) % ny
              nz = Front % elem(e) % nz

              ! Distance from c1 to intersection
              dsc1 = (  (Front % elem(e) % xe - grid % xc(c1)) * nx     &
                      + (Front % elem(e) % ye - grid % yc(c1)) * ny     &
                      + (Front % elem(e) % ze - grid % zc(c1)) * nz  )  &
                   / (lx * nx + ly * ny + lz * nz)

              ! Intersection point
              xs = grid % xc(c1) + dsc1 * lx
              ys = grid % yc(c1) + dsc1 * ly
              zs = grid % zc(c1) + dsc1 * lz

              ! Check if the intersection point is at the element
              do i_ver = 1, Front % elem(e) % nv
                j_ver = i_ver + 1
                if(j_ver > Front % elem(e) % nv) j_ver = 1
                i = Front % elem(e) % v(i_ver)
                j = Front % elem(e) % v(j_ver)
                vec_i(1) = Front % Vert(i) % x_n - xs
                vec_i(2) = Front % Vert(i) % y_n - ys
                vec_i(3) = Front % Vert(i) % z_n - zs
                vec_j(1) = Front % Vert(j) % x_n - xs
                vec_j(2) = Front % Vert(j) % y_n - ys
                vec_j(3) = Front % Vert(j) % z_n - zs
                vec_ixj = Math_Mod_Cross_Product(vec_i, vec_j)
                if( dot_product(vec_ixj(1:3), (/lx,ly,lz/)) < 0.0 ) goto 1
              end do  ! i_ver

              n_int = n_int + 1
              grid % xs(s) = xs
              grid % ys(s) = ys
              grid % zs(s) = zs

  1           continue

              n_fac = n_fac + 1
              Front % face_at_elem(i_cel,s) = e

              ! PRINT '(A,3F12.3)', 'SURFACE FOUND AT C1',  &
              !       GRID % XS(S), GRID % YS(S), GRID % ZS(S)
            end if  ! something with angle
          end if    ! e .ne. 0
        end do      ! i_cel

!       if(n_int < 1) then
!         print *, '# Very strange, no intersections found at face', s
!         print *, '# phi % n(c1) = ', phi % n(c1)
!         print *, '# phi % n(c2) = ', phi % n(c2)
!         print *, '# Front % cell_at_elem(c1) = ', Front % cell_at_elem(c1)
!         print *, '# Front % cell_at_elem(c2) = ', Front % cell_at_elem(c2)
!         print *, '# This error is critical, stopping!'
!         stop
!       end if

        if(n_int > 1) then
          print *, '# Very strange, multiple intersections found at face', s
          print *, '# phi % n(c1) = ', phi % n(c1)
          print *, '# phi % n(c2) = ', phi % n(c2)
          print *, '# This error is critical, stopping!'
          stop
        end if

      end if  ! face crosses 0.5
    end do    ! through faces

    if(Front % n_elems .ne. n_fac) then
      print *, '# It seems that not all face intersections have been found!'
      print *, '# Front % n_elems = ', Front % n_elems
      print *, '# n_fac           = ', n_fac
      print *, '# This error is critical, stopping!'
      stop
    end if

  end if      ! mass transfer

  end subroutine
