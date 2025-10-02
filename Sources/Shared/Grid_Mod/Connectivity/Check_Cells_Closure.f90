!==============================================================================!
  subroutine Check_Cells_Closure(Grid)
!------------------------------------------------------------------------------!
!>  Check integrity of cells' enclosures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, i_fac, s, sh, m, n
  real                 :: dx, dy, dz, sx, sy, sz
  integer, allocatable :: fail_count(:)
!==============================================================================!

  call Profiler % Start('Check_Cells_Closure')

  !------------------------------------------!
  !                                          !
  !   Check integrity of cells' enclosures   !
  !                                          !
  !------------------------------------------!
  allocate(fail_count(-Grid % n_bnd_cells:Grid % n_cells))
  fail_count(:) = 0

  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    Assert(Grid % cells_n_faces(c) .ne. 0)
    do i_fac = 1, Grid % cells_n_faces(c)

      !-----------------------------------!
      !   Take the face and its metrics   !
      !-----------------------------------!
      s = Grid % cells_f(i_fac, c)
      Assert(s .ne. 0)               ! this is just for kicks
      Assert(s .le. Grid % n_faces)  ! this is to remind you that there ...
                                     ! ... are no shadow faces in cells_f
      m = Grid % cells_n_nodes(c)
      n = Grid % faces_n_nodes(s)

      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      sh = Grid % faces_s(s)         ! shadow

      !-------------------------------------!
      !   Check for enclosure with shadow   !
      !-------------------------------------!

      ! If face's nodes are not among cell's nodes, it should be a shadow
      fail_count = 0
      if(.not. Grid % Is_Face_In_Cell(s, c)) then

        Assert(sh .ne. 0)
        Assert(Grid % Is_Face_In_Cell(sh, c))

        ! Surface vector
        call Grid % Faces_Surface(sh, sx, sy, sz)

        ! Ray from cell center to shadow face
        dx = Grid % xf(sh) - Grid % xc(c)
        dy = Grid % yf(sh) - Grid % yc(c)
        dz = Grid % zf(sh) - Grid % zc(c)

        ! This test is stupid, carelessly coppied from case without shadows
        ! because both c1 and c2 are on the same side from the shadow face
        ! if(c .eq. c1) then                  ! (sur)face sign is good
        !   if(dx*sx + dy*sy + dz*sz < 0.0) fail_count = fail_count + 1
        ! else if(c .eq. c2) then             ! (sur)face sign is inverted
        !   if(dx*sx + dy*sy + dz*sz > 0.0) fail_count = fail_count + 1
        ! end if
        if(dx*sx + dy*sy + dz*sz < 0.0) fail_count(c) = fail_count(c) + 1
      end if

      !----------------------------------------!
      !   Check for enclosure without shadow   !
      !----------------------------------------!

      ! If face's nodes are among cell's nodes, it shouldn't be a shadow
      fail_count = 0
      if(Grid % Is_Face_In_Cell(s, c)) then

        ! Surface vector
        call Grid % Faces_Surface(s, sx, sy, sz)

        ! Ray from cell center to face
        dx = Grid % xf(s) - Grid % xc(c)
        dy = Grid % yf(s) - Grid % yc(c)
        dz = Grid % zf(s) - Grid % zc(c)

        if(c .eq. c1) then                  ! (sur)face sign is good
          if(dx*sx + dy*sy + dz*sz < 0.0) fail_count(c) = fail_count(c) + 1
        else if(c .eq. c2) then             ! (sur)face sign is inverted
          if(dx*sx + dy*sy + dz*sz > 0.0) fail_count(c) = fail_count(c) + 1
        end if
      end if
      if(fail_count(c) .gt. 0) call Message % Error(60,                       &
                                   'Cell enclosure without shadows failed!',  &
                                    file=__FILE__)

    end do  ! through cells' faces
  end do    ! through cells
  if(maxval(fail_count) .gt. 0) then
    call Grid % Save_Debug_Vtu("fail",                        &
                               scalar_cell=real(fail_count),  &
                               scalar_name="fail")
    call Message % Error(60,                       &
                        'Cell enclosure with faces test failed!',  &
                         file=__FILE__)
  end if

  call Profiler % Stop('Check_Cells_Closure')

  end subroutine
