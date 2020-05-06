!==============================================================================!
  subroutine Interface_Mod_Create(inter, grid, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)    :: inter(MD, MD)
  type(Grid_Type), target :: grid(MD)
  integer                 :: n_dom

  integer           :: d_1, d_2, k, nks, s, c2, n_1, n_2, s_1, s_2
  real              :: x1, y1, z1, x2, y2, z2, dist, min_dist
  logical           :: found
  character(len=80) :: keys(128)
!==============================================================================!

  !--------------------------------------------!
  !                                            !
  !   Browse through all domain combinations   !
  !                                            !
  !--------------------------------------------!
  do d_1 = 1, n_dom
    do d_2 = 1, n_dom
      inter(d_1, d_2) % n_faces = 0
      if(d_1 .ne. d_2) then

        !-------------------------------------------------------------------!
        !   Try to find specified interface condition between d_1 and d_2   !
        !-------------------------------------------------------------------!
        call Control_Mod_Position_At_Three_Keys('INTERFACE_CONDITION',    &
                                                trim(problem_name(d_1)),  &
                                                trim(problem_name(d_2)),  &
                                                found,                    &
                                                verbose=.false.)

        !-----------------------------------------------------!
        !   Found specification between domains d_1 and d_2   !
        !-----------------------------------------------------!
        if(found) then
          call Control_Mod_Read_Strings_On('BOUNDARY_CONDITIONS',  &
                                           keys, nks, .false.)
          do k = 1, nks
            call To_Upper_Case(keys(k))
          end do

          !---------------------------------------!
          !   Count boundary cells on each side   !
          !---------------------------------------!
          n_1 = 0
          do s = 1, grid(d_1) % n_faces
            c2 = grid(d_1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d_1), c2) .eq. keys(1)) then
                n_1 = n_1 + 1
              end if
            end if
          end do

          n_2 = 0
          do s = 1, grid(d_2) % n_faces
            c2 = grid(d_2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d_2), c2) .eq. keys(2)) then
                n_2 = n_2 + 1
              end if
            end if
          end do

          if(n_1 .ne. n_2) then
            print *, '# Number of cells at the interface between ',  &
                      trim(problem_name(d_1)), ' and ',               &
                      trim(problem_name(d_2)), ' is not the same!'
            print *, '# Only conformal mappings are supported.  Exiting!'
            stop
          else
            print '(a,a,a,a,a,i6,a)', ' # Domains ', trim(problem_name(d_1)),  &
                                            ' and ', trim(problem_name(d_2)),  &
                      ' are connected with ', n_1, ' interface cells!'
          end if

          !---------------------------------------!
          !   Store boundary cells on each side   !
          !---------------------------------------!
          inter(d_1, d_2) % n_faces = n_1  ! the same as n_2
          allocate(inter(d_1, d_2) % faces_1(n_1))
          allocate(inter(d_1, d_2) % faces_2(n_2))
          allocate(inter(d_1, d_2) % close_in_1(n_1))
          allocate(inter(d_1, d_2) % close_in_2(n_2))
          inter(d_1, d_2) % faces_1(:) = 0
          inter(d_1, d_2) % faces_2(:) = 0
          inter(d_1, d_2) % close_in_1(:) = 0
          inter(d_1, d_2) % close_in_2(:) = 0

          ! On the side of d_1
          n_1 = 0
          do s = 1, grid(d_1) % n_faces
            c2 = grid(d_1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d_1), c2) .eq. keys(1)) then
                n_1 = n_1 + 1
                inter(d_1, d_2) % faces_1(n_1) = s
              end if
            end if
          end do

          ! On the side of d_2
          n_2 = 0
          do s = 1, grid(d_2) % n_faces
            c2 = grid(d_2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d_2), c2) .eq. keys(2)) then
                n_2 = n_2 + 1
                inter(d_1, d_2) % faces_2(n_2) = s
              end if
            end if
          end do
        end if
      end if
    end do
  end do

  !-----------------------------------------------------!
  !   Find closest faces from the neighbouring domain   !
  !-----------------------------------------------------!
  do d_1 = 1, n_dom
    do d_2 = 1, n_dom
      if(inter(d_1, d_2) % n_faces > 0) then

        do n_1 = 1, inter(d_1, d_2) % n_faces
          s_1 = inter(d_1, d_2) % faces_1(n_1)
          x1 = grid(d_1) % xf(s_1)
          y1 = grid(d_1) % yf(s_1)
          z1 = grid(d_1) % zf(s_1)
          min_dist = HUGE
          do n_2 = 1, inter(d_1, d_2) % n_faces
            s_2 = inter(d_1, d_2) % faces_2(n_2)
            x2 = grid(d_2) % xf(s_2)
            y2 = grid(d_2) % yf(s_2)
            z2 = grid(d_2) % zf(s_2)

            dist = Math_Mod_Distance_Squared(x1, y1, z1, x2, y2, z2)
            if(dist < min_dist) then
              inter(d_1, d_2) % close_in_2(n_1) = n_2
              inter(d_1, d_2) % close_in_1(n_2) = n_1
              min_dist = dist
            end if

          end do
        end do
      end if
    end do
  end do

! ! Write some debugging info
! do d_1 = 1, n_dom
!   do d_2 = 1, n_dom
!     if(inter(d_1, d_2) % n_faces > 0) then
!       do n_1 = 1, inter(d_1, d_2) % n_faces
!         s_1 = inter(d_1, d_2) % faces_1(n_1)
!         n_2 = inter(d_1, d_2) % close_in_2(n_1)
!         s_2 = inter(d_1, d_2) % faces_2(n_2)
!         WRITE(100,'(6F9.4)')  &
!           grid(d_1) % xf(s_1), grid(d_1) % yf(s_1), grid(d_1) % zf(s_1),  &
!           grid(d_2) % xf(s_2), grid(d_2) % yf(s_2), grid(d_2) % zf(s_2)
!       end do
!     end if
!   end do
! end do

  end subroutine
