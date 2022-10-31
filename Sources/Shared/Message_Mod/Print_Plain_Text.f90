!==============================================================================!
  subroutine Print_Plain_Text(Msg, width, message_text)
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)      :: Msg
  integer,      intent(in) :: width
  character(*), intent(in) :: message_text
!-----------------------------------[Locals]-----------------------------------!
  type(Tokenizer_Type) :: Tok
  integer              :: i
  integer              :: cur_p, nex_p
  character(DL)        :: line
!==============================================================================!

  !---------------------------------------------------!
  !                                                   !
  !   Load the argument into tokenizer and parse it   !
  !                                                   !
  !---------------------------------------------------!
  Tok % whole = message_text
  call Tok % Parse()

  ! Initialize the line
  line      = ' '
  line(2:2) = '#'
  cur_p = 4

  !-----------------------------------------------------!
  !                                                     !
  !   Browse through tokens and print them one by one   !
  !                                                     !
  !-----------------------------------------------------!
  i = 1
  do while (i .le. Tok % n_tokens)
1   continue
    nex_p = cur_p + len_trim(Tok % tokens(i))

    !-------------------------------------------------!
    !   Keep on filling up this line, it still fits   !
    !-------------------------------------------------!
    if(nex_p < width .and. Tok % tokens(i) .ne. '\n') then
      write(line(cur_p:nex_p), '(a)')  trim(Tok % tokens(i)) // ' '

    !---------------------------------------------------------!
    !   Line is too long, print it and reset for new tokens   !
    !---------------------------------------------------------!
    else

      ! Print what you have up to now
      print '(a)', trim(line)

      ! (Re)initialize the line
      line      = ' '
      line(2:2) = '#'
      cur_p = 4

      ! If you came here because of new line character, skip it.
      if(Tok % tokens(i) .eq. '\n') i = i + 1
      goto 1
    end if

    cur_p = nex_p + 1
    i = i + 1
  end do

  !--------------------------------------------------------!
  !                                                        !
  !   Print the last line which was still in the forming   !
  !                                                        !
  !--------------------------------------------------------!
  print '(a)', trim(line)

  end subroutine
