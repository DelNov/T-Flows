!==============================================================================!
  subroutine Parse(Tok)
!------------------------------------------------------------------------------!
!>  Breaks down a string contained in an instance of Tokenizer_Type into
!>  separate tokens. It is useful for text processing where extracting
!>  individual words or symbols from a line is necessary, such as when
!>  importing grids in various formats, reading control file and alike.
!------------------------------------------------------------------------------!
!   Functionality:
!
!   * Initial Setup:
!     - Determines the length of the string to be tokenized and stores the
!       first and last characters of the string in Tok % first and Tok % last,
!       respectively.
!   * Tokenization Logic:
!     - Initializes the token count (n_tokens) to zero.
!     - Identifies tokens by scanning through the string: A token starts when
!       a character is a non-space ('!' or higher in ASCII) following a space
!       or at the beginning of the string. The end of a token is marked by a
!       space following a non-space character.
!     - Updates Tok % s and Tok % e arrays with the start and end positions
!       of each token.
!     - Uses assertions to ensure the number of tokens does not exceed the
!       predefined maximum (MAX_TOKENS).
!   * Extracting Tokens:
!     - Extracts individual tokens from the whole string based on identified
!       start and end positions and stores them in the tokens array of Tok.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Tokenizer_Type), intent(inout) :: Tok  !! parent class containing the
                                               !! string to be tokenized
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, l
!==============================================================================!

  ! Fetch the first and the last character
  l = len_trim(Tok % whole)
  Tok % first = Tok % whole(1:1)
  Tok % last  = Tok % whole(l:l)

  ! Then parse
  Tok % n_tokens = 0
  if(Tok % whole(1:1) >= '!') then
    Tok % n_tokens = 1
    Tok % s(1)=1
  end if
  do i = 1, len(Tok % whole) - 1
    if( Tok % whole(i:  i  ) <  '!' .and.  &
        Tok % whole(i+1:i+1) >= '!') then
      Tok % n_tokens = Tok % n_tokens + 1
      Assert(Tok % n_tokens <= MAX_TOKENS)
      Tok % s(Tok % n_tokens) = i+1
    end if
    if( Tok % whole(i  :i  ) >= '!' .and.  &
        Tok % whole(i+1:i+1) <  '!') then
      Tok % e(Tok % n_tokens) = i
    end if
  end do

  ! Chop them up
  do i = 1, Tok % n_tokens
    Tok % tokens(i) = Tok % whole(Tok % s(i) : Tok % e(i))
  end do

  end subroutine
