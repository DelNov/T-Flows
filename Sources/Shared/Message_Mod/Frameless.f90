!==============================================================================!
  subroutine Frameless(Msg, width, message_text)
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)      :: Msg
  integer,      intent(in) :: width
  character(*), intent(in) :: message_text
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i
  integer       :: cur_p, nex_p
  character(DL) :: out_line
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Msg)
!==============================================================================!

  !---------------------------------------------------!
  !                                                   !
  !   Load the argument into tokenizer and parse it   !
  !                                                   !
  !---------------------------------------------------!
  Line % whole = message_text
  call Line % Parse()

  ! Initialize the output line
  out_line      = ' '
  out_line(2:2) = '#'
  cur_p         = 4    ! starts at 4 to allow for one leading space in output

  !-----------------------------------------------------!
  !                                                     !
  !   Browse through tokens and print them one by one   !
  !                                                     !
  !-----------------------------------------------------!
  i = 1
  do while (i .le. Line % n_tokens)
1   continue
    nex_p = cur_p + len_trim(Line % tokens(i))

    !-------------------------------------------------!
    !   Keep on filling up this line, it still fits   !
    !   (+3 here is to compensate you started at 4)   !
    !-------------------------------------------------!
    if(nex_p < width + 3 .and. Line % tokens(i) .ne. '\n') then
      write(out_line(cur_p:nex_p), '(a)')  trim(Line % tokens(i)) // ' '

    !---------------------------------------------------------!
    !   Line is too long, print it and reset for new tokens   !
    !---------------------------------------------------------!
    else

      ! Print what you have up to now
      print '(a)', trim(out_line)

      ! (Re)initialize the output line
      out_line      = ' '
      out_line(2:2) = '#'
      cur_p         = 4

      ! If you came here because of new line character, skip it.
      if(Line % tokens(i) .eq. '\n') i = i + 1
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
  print '(a)', trim(out_line)

  end subroutine
