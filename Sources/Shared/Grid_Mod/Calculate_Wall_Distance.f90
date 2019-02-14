!==============================================================================!
  subroutine Grid_Mod_Calculate_Wall_Distance(grid)
!------------------------------------------------------------------------------!
!   Calculate distance from the cell center to the nearest wall.               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include '../Shared/Approx_Real.int'
!----------------------------------[Calling]-----------------------------------!
  real :: Distance_Squared
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c1, c2
  integer              :: n_wall_colors
  integer, allocatable :: wall_colors(:)
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
  print *, '# This is needed for RANS and HYBRID turbulence models.            '
  print *, '#------------------------------------------------------------------'
  call Tokenizer_Mod_Read_Line(5)
  n_wall_colors = line % n_tokens
  allocate(wall_colors(n_wall_colors))
  do b = 1, n_wall_colors
    read(line % tokens(b), *) wall_colors(b)
  end do

  if( (n_wall_colors .eq. 1) .and. (wall_colors(1) .eq. 0) ) then
    grid % wall_dist = 1.0
    print *, '# Distance to the wall set to 1.0 everywhere !'
  else
    do c1 = 1, grid % n_cells
      if(mod(c1,10000) .eq. 0) then
        write(*,'(a2, f5.0, a14)') ' #', (100.*c1/(1.*grid % n_cells)),  &
                                   ' % complete...'
      end if
      do b = 1, n_wall_colors
        do c2 = grid % bnd_cond % color_f( wall_colors(b) ),  &
                grid % bnd_cond % color_l( wall_colors(b) ),  &
                -1
          grid % wall_dist(c1) = min(grid % wall_dist(c1),            &
                                     Distance_Squared(grid % xc(c1),  &
                                                      grid % yc(c1),  &
                                                      grid % zc(c1),  &
                                                      grid % xc(c2),  &
                                                      grid % yc(c2),  &
                                                      grid % zc(c2)))
        end do
      end do
    end do

    grid % wall_dist(:) = sqrt(grid % wall_dist(:))

    print *, '# Distance to the wall calculated !'
  end if

  deallocate(wall_colors)

  end subroutine
