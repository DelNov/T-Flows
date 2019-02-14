!==============================================================================!
  subroutine Connect_Domains(grid)
!------------------------------------------------------------------------------!
!   Connects two problem domains, one with periodic streamwise boundary        !
!   conditions and another with inflow-outflow.                                !
!                                                                              !
!   Note:                                                                      !
!                                                                              !
!   Situations like the on depicted bellow are now working.                    !
!                                                                              !
!   +-------+   +_                                                             !
!   |       |   | ~-_                                                          !
!   |   o---|   |    ~-_                                                       !
!   |       |   |---o   +                                                      !
!   +-------+   +_      |                                                      !
!                 ~-_   |                                                      !
!                    ~-_|                                                      !
!                       +                                                      !
!                                                                              !
!   Constrains:                                                                !
!                                                                              !
!   First domain must be the channel-like, periodic in streamwise direction.   !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,     only: ONE_THIRD, TWO_THIRDS
  use Tokenizer_Mod, only: Tokenizer_Mod_Read_Line, line
  use Grid_Mod,      only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  include '../Shared/Approx_Real.int'
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, c1, c11, c12, c21, c22, s1, s2
  integer           :: color_copy,  x_copy, y_copy, z_copy
  real              :: xc_12, xc_22
  real              :: yc_12, yc_22
  real              :: zc_12, zc_22
  real, parameter   :: SMALL = 1.0e-4
  character(len=80) :: answer
!==============================================================================!

  grid % n_copy = 0
  x_copy = 0
  y_copy = 0
  z_copy = 0

1 print *, '#======================================================'
  print *, '# Type ordinal number(s) of inflow boundary(s) that    '
  print *, '# will be linked to precursor domain(s) (skip to exit):'
  print *, '#------------------------------------------------------'
  call Tokenizer_Mod_Read_Line(5)
  read(line % tokens(1), *) answer
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') then
    return
  else 
    read(line % tokens(1), *) color_copy 
  end if    

  !-------!
  !   X   !
  !-------! 
  do s1 = 1, grid % n_faces
    if(mod(s1,100000) .eq. 0) then
      print *, ((ONE_THIRD*s1)/(1.0*grid % n_faces)) * 100, 'Complete'
    end if
    c11 = grid % faces_c(1,s1)
    c12 = grid % faces_c(2,s1)
    if( abs(grid % dx(s1)) > 1.0e-9 ) then
      do s2 = 1, grid % n_faces 
        c21 = grid % faces_c(1,s2)
        c22 = grid % faces_c(2,s2)
        if(c22 < 0) then
          if(grid % bnd_cond % color(c22) .eq. color_copy) then

            yc_12 = 0.0
            zc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              yc_12 = yc_12 + grid % yn(grid % faces_n(i,s1))
              zc_12 = zc_12 + grid % zn(grid % faces_n(i,s1))
            end do
            yc_12 = yc_12 / (real(grid % faces_n_nodes(s1)))
            zc_12 = zc_12 / (real(grid % faces_n_nodes(s1)))

            yc_22 = 0.0
            zc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              yc_22 = yc_22 + grid % yn(grid % faces_n(i,s2))
              zc_22 = zc_22 + grid % zn(grid % faces_n(i,s2))
            end do
            yc_22 = yc_22 / (real(grid % faces_n_nodes(s2)))
            zc_22 = zc_22 / (real(grid % faces_n_nodes(s2)))

            if( Approx_Real( yc_22, yc_12, tol=SMALL ) .and. &
                Approx_Real( zc_22, zc_12, tol=SMALL ) ) then
              grid % n_copy = grid % n_copy + 1
              x_copy = x_copy + 1
              if( abs(grid % xc(c11)-grid % xc(c22)) <  &
                  abs(grid % xc(c12)-grid % xc(c22))) c1 = c11
              if( abs(grid % xc(c11)-grid % xc(c22)) >  &
                  abs(grid % xc(c12)-grid % xc(c22))) c1 = c12
              grid % bnd_cond % copy_s(1, grid % n_copy) = c1
              grid % bnd_cond % copy_s(2, grid % n_copy) = c21  ! inside domain
              grid % bnd_cond % copy_c(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  !-------!
  !   Y   !
  !-------! 
  do s1 = 1, grid % n_faces
    if(mod(s1,100000) .eq. 0) then
      print *, (ONE_THIRD + (ONE_THIRD*s1)  &
               / (1.0 * grid % n_faces)) * 100.0, 'Complete'
    end if
    c11 = grid % faces_c(1,s1)
    c12 = grid % faces_c(2,s1)
    if( abs(grid % dy(s1)) > 1.0e-9 ) then
      do s2 = 1, grid % n_faces 
        c21 = grid % faces_c(1,s2)
        c22 = grid % faces_c(2,s2)
        if(c22 < 0) then
          if(grid % bnd_cond % color(c22) .eq. color_copy) then

            xc_12 = 0.0
            zc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              xc_12 = xc_12 + grid % xn(grid % faces_n(i,s1))
              zc_12 = zc_12 + grid % zn(grid % faces_n(i,s1))
            end do
            xc_12 = xc_12 / (real(grid % faces_n_nodes(s1)))
            zc_12 = zc_12 / (real(grid % faces_n_nodes(s1)))

            xc_22 = 0.0
            zc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              xc_22 = xc_22 + grid % xn(grid % faces_n(i,s2))
              zc_22 = zc_22 + grid % zn(grid % faces_n(i,s2))
            end do
            xc_22 = xc_22 / (real(grid % faces_n_nodes(s2)))
            zc_22 = zc_22 / (real(grid % faces_n_nodes(s2)))

            if( Approx_Real( xc_22, xc_12, tol=SMALL ) .and. &
                Approx_Real( zc_22, zc_12, tol=SMALL ) ) then
              grid % n_copy = grid % n_copy + 1 
              y_copy = y_copy + 1
              if( abs(grid % yc(c11)-grid % yc(c22)) <  &
                  abs(grid % yc(c12)-grid % yc(c22))) c1 = c11
              if( abs(grid % yc(c11)-grid % yc(c22)) >  &
                  abs(grid % yc(c12)-grid % yc(c22))) c1 = c12
              grid % bnd_cond % copy_s(1, grid % n_copy) = c1
              grid % bnd_cond % copy_s(2, grid % n_copy) = c21  ! inside domain
              grid % bnd_cond % copy_c(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  !-------!
  !   Z   !
  !-------! 
  do s1 = 1, grid % n_faces
    if(mod(s1,100000) .eq. 0) then
      print *, (TWO_THIRDS + (ONE_THIRD*s1)  &
               / (1.0*grid % n_faces)) * 100.0, 'Complete'
    end if
    c11 = grid % faces_c(1,s1)
    c12 = grid % faces_c(2,s1)
    if( abs(grid % dz(s1)) > 1.0e-9 ) then
      do s2 = 1, grid % n_faces 
        c21 = grid % faces_c(1,s2)
        c22 = grid % faces_c(2,s2)
        if(c22 < 0) then
          if(grid % bnd_cond % color(c22) .eq. color_copy) then

            yc_12 = 0.0
            xc_12 = 0.0
            do i=1,grid % faces_n_nodes(s1)
              yc_12 = yc_12 + grid % yn(grid % faces_n(i,s1))
              xc_12 = xc_12 + grid % xn(grid % faces_n(i,s1))
            end do
            yc_12 = yc_12 / (real(grid % faces_n_nodes(s1)))
            xc_12 = xc_12 / (real(grid % faces_n_nodes(s1)))

            yc_22 = 0.0
            xc_22 = 0.0
            do i=1,grid % faces_n_nodes(s2)
              yc_22 = yc_22 + grid % yn(grid % faces_n(i,s2))
              xc_22 = xc_22 + grid % xn(grid % faces_n(i,s2))
            end do
            yc_22 = yc_22 / (real(grid % faces_n_nodes(s2)))
            xc_22 = xc_22 / (real(grid % faces_n_nodes(s2)))

            if( Approx_Real( yc_22, yc_12, tol=SMALL ) .and. &
                Approx_Real( xc_22, xc_12, tol=SMALL ) ) then
              grid % n_copy = grid % n_copy + 1 
              z_copy = z_copy + 1
              if( abs(grid % zc(c11)-grid % zc(c22)) <  &
                  abs(grid % zc(c12)-grid % zc(c22))) c1 = c11
              if( abs(grid % zc(c11)-grid % zc(c22)) >  &
                  abs(grid % zc(c12)-grid % zc(c22))) c1 = c12
              grid % bnd_cond % copy_s(1, grid % n_copy) = c1
              grid % bnd_cond % copy_s(2, grid % n_copy) = c21  
              grid % bnd_cond % copy_c(c22) = c1
            end if
          end if
        end if
      end do
    end if
  end do

  print *, '# n copy = ', grid % n_copy
  print *, '# x copy = ', x_copy
  print *, '# x copy = ', y_copy
  print *, '# x copy = ', z_copy

  goto 1

  end subroutine
