!==============================================================================!
  subroutine Particles_Mod_Find_Nearest_Cell(flow, part)
!------------------------------------------------------------------------------!
!   Finds a cell closest to a particle.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod, only: Field_Type
  use Grid_Mod,  only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Particle_Type)         :: part
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2, c
  real                     :: xf, yf, zf  ! face center coordinates
  real                     :: xc, yc, zc  ! cell center coordinates
  real                     :: d_sq        ! distance squared
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  !----------------------------------------------!
  !                                              !
  !   First level comments are framed with one   !
  !   extra empty row above and below text       !
  !                                              !
  !----------------------------------------------!

  !------------------------------------------------------!
  !                                                      ! 
  !   First level comments are written like sentences,   ! 
  !   first word capital, and can spread over several    !
  !   lines like this very commet.                       !
  !                                                      !
  !------------------------------------------------------!

  !---------------------------!
  !   Second level comments   !
  !---------------------------!

  !-----------------------------------------------------!
  !     First and second level comments have frames     !
  !   three spaces thick to the left and to the right   !
  !-----------------------------------------------------!
                                                    !123! <- see thickness

  ! Third level comments follow sentence rules, starting with upper-case letter

  ! If a third level comment is too long ...
  ! ... it can spread over several lines ...
  ! ... with dots for continuation.

  !--------------------------------------------!
  !   Browse through all faces as an example   !
  !--------------------------------------------!
  do s = 1, grid % n_faces

    ! Take indices of cells surrounding the face
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(1, s)

    ! Take cell coordinates
    xf = grid % xf(s)
    yf = grid % yf(s)
    zf = grid % zf(s)

  end do

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  do c = 1, grid % n_cells

    ! Take cell centre
    xc = grid % xc(c)
    yc = grid % yc(c)
    zc = grid % zc(c)

    ! Distance squared from the particle to cell centre
    d_sq = (xc - part % x)**2  &
         + (yc - part % y)**2  &
         + (zc - part % z)**2


  end do

  end subroutine
