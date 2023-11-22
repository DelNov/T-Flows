!==============================================================================!
  subroutine Parse(Tok)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Tokenizer_Type), intent(inout) :: Tok
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
