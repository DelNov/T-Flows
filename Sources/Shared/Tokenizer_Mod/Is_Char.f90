
!==============================================================================!
  pure logical function Is_Char(Tok, ch)
!------------------------------------------------------------------------------!
!>  Checks whether a character should be treated as part of a token when
!>  tokenizing strings.
!>
!>  Normally, printable characters are considered token characters, while
!>  control characters and blanks are treated as separators.  However, the ANSI
!>  escape character ESC is treated as a token character as well.  This allows
!>  coloured terminal strings, such as ESC//"[96m", to pass through the
!>  tokenizer without losing the ESC character.
!>
!>  Note that ANSI colour sequences are still counted as ordinary characters
!>  when estimating string length.  This slightly over-estimates the visible
!>  length of coloured strings, but this is acceptable for formatted diagnostic
!>  messages.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Tokenizer_Type), intent(in)   :: Tok  !! parent class containing the
                                              !! string to be tokenized
  character(len=1), intent(in)        :: ch   !! input character
!==============================================================================!

  ! Normal printable characters are token characters.
  ! ANSI ESC is also a token character, even though ASCII 27 < '!'.
  Is_Char = iachar(ch) >= iachar('!') .or. iachar(ch) == 27

  end function

