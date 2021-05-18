!==============================================================================!
  subroutine Find_Periodic_Faces(Grid)
!------------------------------------------------------------------------------!
!   Periodic faces are needed for Lagrangian particle tracking.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, run    ! counters
!==============================================================================!

  !----------------------------------------!
  !   Count the number of periodic faces   !
  !----------------------------------------!
  do run = 1, 2

    ! Initialize number of periodic faces
    if(run .eq. 1) Grid % n_per_faces = 0
    if(run .eq. 2) then
      allocate(Grid % per_faces(Grid % n_per_faces))
      Grid % n_per_faces = 0
    end if

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      ! In x direction
      if(c2 .gt. 0) then

        ! In x direction
        if( (abs(Grid % per_x) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(Grid % dx(s))                    &
                                  + abs(Grid % xc(c1) - Grid % xc(c2)),  &
                                    Grid % per_x) ) then
          Grid % n_per_faces = Grid % n_per_faces + 1
          if(run .eq. 2) then
            Grid % per_faces(Grid % n_per_faces) = s  ! store periodic face
            if(Grid % xc(c2) > Grid % xc(c1))                      &
              Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 1
          end if
        end if

        ! In y direction
        if( (abs(Grid % per_y) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(Grid % dy(s))                    &
                                  + abs(Grid % yc(c1) - Grid % yc(c2)),  &
                                    Grid % per_y) ) then
          Grid % n_per_faces = Grid % n_per_faces + 1
          if(run .eq. 2) then
            Grid % per_faces(Grid % n_per_faces) = s
            if(Grid % yc(c2) > Grid % yc(c1))                      &
              Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 2
          end if
        end if

        ! In z direction
        if( (abs(Grid % per_z) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(Grid % dz(s))                    &
                                  + abs(Grid % zc(c1) - Grid % zc(c2)),  &
                                    Grid % per_z) ) then
          Grid % n_per_faces = Grid % n_per_faces + 1
          if(run .eq. 2) then
            Grid % per_faces(Grid % n_per_faces) = s
            if(Grid % zc(c2) > Grid % zc(c1))                      &
              Grid % bnd_cond % color(s) = Grid % n_bnd_cond + 3
          end if
        end if
      end if
    end do
  end do

  end subroutine
