!==============================================================================!
!                                  hexahedron                                  !
!                                                                              !
!                      7-----------8             7-----------8                 !
!                     /|          /|            /|          /|                 !
!                    /           / |           /    (6)    / |                 !
!                   /  |    (4) /  |          /  |        /  |                 !
!                  5-----------6   |         5-----------6   |                 !
!                  |(5)|       |(3)|         |   |       |   |                 !
!                  |   3- - - -|- -4         |   3- - - -|- -4                 !
!                  |  / (2)    |  /          |  /        |  /                  !
!                  |           | /           |      (1)  | /                   !
!                  |/          |/            |/          |/                    !
!                  1-----------2             1-----------2                     !
!                                                                              !
!                  (3) and (4) are behind                                      !
!------------------------------------------------------------------------------!
!   This is the native block and cell numbering for Generator.                 !
!                                                                              !
!   Unfortunatelly, it is not the same as the native Gambit's neutral file     !
!   format, which is used in Convertor - the face order is different.          !
!                                                                              !
!   Changing it in order to make it consistent is not very straighforward      !
!   since blocks are connected using this mapping too, meaning that if it      !
!   changes, all domain files should change accordingly. Doesn't make sense.   !
!==============================================================================!

  ! Declarations
  integer, dimension(6,4) :: hex_block

  hex_block = transpose(reshape( (/ 1, 2, 4, 3,  &
                                    1, 5, 6, 2,  &
                                    2, 6, 8, 4,  &
                                    4, 8, 7, 3,  &
                                    3, 7, 5, 1,  &
                                    5, 7, 8, 6  /), (/4, 6/) ))
