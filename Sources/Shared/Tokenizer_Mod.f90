!==============================================================================!
  module Tokenizer_Mod
!------------------------------------------------------------------------------!
!   Holds functionality to parse lines into tokens.  Tokenizes lines.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_TOKENS = QL / 2
!==============================================================================!

  !--------------------!
  !   Tokenizer type   !
  !--------------------!
  type Tokenizer_Type
    character(QL) :: whole               ! whole string
    character(SL) :: tokens(MAX_TOKENS)  ! tokens
    integer       :: n_tokens            ! number of tokens
    integer       :: s(MAX_TOKENS),  &   ! tokens starts ...
                     e(MAX_TOKENS)       ! ... and ends
    character(1)  :: first, last         ! first and last characters in whole

    contains
      procedure :: Parse

  end type

  contains

  include 'Tokenizer_Mod/Parse.f90'

  end module
