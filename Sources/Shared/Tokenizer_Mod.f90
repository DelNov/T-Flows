#include "../Shared/Assert.h90"

!==============================================================================!
  module Tokenizer_Mod
!------------------------------------------------------------------------------!
!   Holds functionality to parse lines into tokens.  Tokenizes lines.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_TOKENS = 2048
!==============================================================================!

  !--------------------!
  !   Tokenizer type   !
  !--------------------!
  type Tokenizer_Type
    character(MAX_TOKENS*2) :: whole               ! whole string
    character(SL)           :: tokens(MAX_TOKENS)  ! tokens
    integer                 :: n_tokens            ! number of tokens
    integer                 :: s(MAX_TOKENS),  &   ! tokens starts ...
                               e(MAX_TOKENS)       ! ... and ends
    character(1)            :: first, last         ! the first and last ...
                                                   ! character in the whole
    contains
      procedure :: Parse

  end type

  contains

#   include "Tokenizer_Mod/Parse.f90"

  end module
