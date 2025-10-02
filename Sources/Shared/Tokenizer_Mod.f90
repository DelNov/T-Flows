#include "../Shared/Assert.h90"

!==============================================================================!
  module Tokenizer_Mod
!------------------------------------------------------------------------------!
!>  Provides functionality to parse lines of text into individual tokens,
!>  facilitating the handling and analysis of textual data.  It is first used
!>  in Message_Mod, then included in File_Mod and, clearly, is available in
!>  all the modules and functions using them.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_TOKENS = 2048  !! maximum number of tokens
                                           !! hat can be handled
!==============================================================================!

  !--------------------!
  !   Tokenizer type   !
  !--------------------!
  !> Includes data members and a procedure for tokenizing strings
  type Tokenizer_Type
    character(MAX_TOKENS*2) :: whole               !! whole string
    character(SL)           :: tokens(MAX_TOKENS)  !! tokens
    integer                 :: n_tokens            !! number of tokens
    integer                 :: s(MAX_TOKENS),  &   !! token starts
                               e(MAX_TOKENS)       !! token ends
    character(1)            :: first, last         !! the first and the last
                                                   !! character in the whole
    contains
      procedure :: Parse  !! procedure to tokenize the string stored in whole

  end type

  !---------------------------!
  !   Singleton object Line   !
  !---------------------------!
  type(Tokenizer_Type) :: Line  !! singleton object Line for the entire T-Flows

  contains

#   include "Tokenizer_Mod/Parse.f90"

  end module
