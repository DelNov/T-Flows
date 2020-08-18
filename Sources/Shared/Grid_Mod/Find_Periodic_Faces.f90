!==============================================================================!
  subroutine Grid_Mod_Find_Periodic_Faces(grid)
!------------------------------------------------------------------------------!
!   Periodic faces are needed for Lagrangian particle tracking.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, run    ! counters
!==============================================================================!

  !----------------------------------------!
  !   Count the number of periodic faces   !
  !----------------------------------------!
  do run = 1, 2

    ! Initialize number of periodic faces
    if(run .eq. 1) grid % n_per_faces = 0
    if(run .eq. 2) then
      allocate(grid % per_faces(grid % n_per_faces))
      grid % n_per_faces = 0
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)

      ! In x direction
      if(c2 .gt. 0) then

        ! In x direction
        if( (abs(grid % per_x) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(grid % dx(s))                    &
                                  + abs(grid % xc(c1) - grid % xc(c2)),  &
                                    grid % per_x) ) then
          grid % n_per_faces = grid % n_per_faces + 1
          if(run .eq. 2) then
            grid % per_faces(grid % n_per_faces) = s  ! store periodic face
            grid % bnd_cond % color(s) = grid % n_bnd_cond + 1
          end if
        end if

        ! In y direction
        if( (abs(grid % per_y) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(grid % dy(s))                    &
                                  + abs(grid % yc(c1) - grid % yc(c2)),  &
                                    grid % per_y) ) then
          grid % n_per_faces = grid % n_per_faces + 1
          if(run .eq. 2) then
            grid % per_faces(grid % n_per_faces) = s
            grid % bnd_cond % color(s) = grid % n_bnd_cond + 2
          end if
        end if

        ! In z direction
        if( (abs(grid % per_z) > NANO) .and.                             &
            Math_Mod_Approx_Real(   abs(grid % dz(s))                    &
                                  + abs(grid % zc(c1) - grid % zc(c2)),  &
                                    grid % per_z) ) then
          grid % n_per_faces = grid % n_per_faces + 1
          if(run .eq. 2) then
            grid % per_faces(grid % n_per_faces) = s
            grid % bnd_cond % color(s) = grid % n_bnd_cond + 3
          end if
        end if
      end if
    end do
  end do

  end subroutine
