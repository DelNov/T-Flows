!==============================================================================!
  subroutine Template_Subroutine(<arguments>)
!------------------------------------------------------------------------------!
!   Template for a subroutine/function.                                        !
!                                                                              !
!   Keep all coding and commenting up to eighty columns wide.  If you go       !
!   beyond that, the code will look terrible when printed, and you do want     !
!   to analyze printed text for better concentration.                          !
!                                                                              !
!   Also note that indentation is by 2!  Subroutine name starts at column 3,   !
!   coding starts in the same column, and next levels of indentation at 5,     !
!   7, 9, etc.                                                                 !
!                                                                              !
!   Note that there is a separator before the function header (this section)   !
!   and line which introduces modules.  If not, the characters "[Modules]"     !
!   get cramped with text written here.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  ...
  ..
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  ...
  ..
!----------------------------------[Calling]-----------------------------------!
  ...
  ..
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: BIG = 1000
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, k
  character(len=SL)    :: name_in
  real, dimension(5,5) :: small_matrix
!------------------------[Avoid unused parent warning]-------------------------!
!==============================================================================!

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

  do s = 1, Grid % n_faces     ! fourth level comments start with small letters

    c1 = Grid % faces_c(1, s)  ! they are too small to qualify as sentences
    c2 = Grid % faces_c(2, s)  ! keep them two columns apart from the last ...
                               ! column with the code
                            !12! <- see distance from the last command
  end do

  end subroutine
