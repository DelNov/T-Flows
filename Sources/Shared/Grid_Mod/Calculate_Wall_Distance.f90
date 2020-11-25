!==============================================================================!
  subroutine Grid_Mod_Calculate_Wall_Distance(grid)
!------------------------------------------------------------------------------!
!   Calculate distance from the cell center to the nearest wall.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c1, c2
  integer              :: n_wall_colors
  integer              :: processed_cells
  integer, allocatable :: wall_colors(:)
  character(SL)        :: answer
!==============================================================================!

  !------------------------------------------------------------------!
  !   Calculate distance from the cell center to the nearest wall.   !
  !------------------------------------------------------------------!
  !     => depends on: xc,yc,zc inside and on the boundary.          !
  !     <= gives:      wall_dist                                     !
  !------------------------------------------------------------------!
  grid % wall_dist = HUGE

  call Grid_Mod_Print_Bnd_Cond_List(grid)
  print *, '#=================================================================='
  print *, '# Calculating distance from the walls                              '
  print *, '#------------------------------------------------------------------'
  print *, '# Type ordinal number(s) of wall or wall_flux boundary condition(s)'
  print *, '# from the boundary condition list (see above) separated by space. '
  print *, '# Cells'' centers distances to the nearest wall will be calculated '
  print *, '# for the listed wall boundary(s).                                 '
  print *, '#                                                                  '
  print *, '# This is needed for RANS and HYBRID turbulence models as well as  '
  print *, '# for proper initialization with potential pressure-like field.    '
  print *, '#                                                                  '
  print *, '# Type skip to skip this and set wall distance to one everywhere.  '
  print *, '#------------------------------------------------------------------'
  call File_Mod_Read_Line(5)
  answer = line % tokens(1)
  call To_Upper_Case(answer)

  ! User wants to skip calculation of wall distance
  if( answer .eq. 'SKIP' ) then
    grid % wall_dist = 1.0
    print *, '# Distance to the wall set to 1.0 everywhere !'

  ! Calculation of wall distance
  else

    n_wall_colors = line % n_tokens
    allocate(wall_colors(n_wall_colors))
    do b = 1, n_wall_colors
      read(line % tokens(b), *) wall_colors(b)
    end do

    processed_cells = 0

    !$omp parallel do
    do c1 = -grid % n_bnd_cells, grid % n_cells

      processed_cells = processed_cells + 1

      ! Write progress and stay in the same line
      ! (achieved with advance='no' and achar(13))
      write(*,'(a2,f5.0,a14,a1)', advance='no')               &
        ' #',                                                 &
        (100. * real(processed_cells)                         &
               / real(grid % n_bnd_cells + grid % n_cells)),  &
        ' % complete...', achar(13)

      do b = 1, n_wall_colors
        do c2 = grid % bnd_cond % color_s_cell( wall_colors(b) ),  &
                grid % bnd_cond % color_e_cell( wall_colors(b) ),  &
                -1
          grid % wall_dist(c1) =  &
            min(grid % wall_dist(c1),                     &
                Math_Mod_Distance_Squared(grid % xc(c1),  &
                                          grid % yc(c1),  &
                                          grid % zc(c1),  &
                                          grid % xc(c2),  &
                                          grid % yc(c2),  &
                                          grid % zc(c2)))
        end do
      end do
    end do
    !$omp end parallel do

    grid % wall_dist(:) = sqrt(grid % wall_dist(:))

    print *, '# Distance to the wall calculated !'

    deallocate(wall_colors)
  end if

  end subroutine
