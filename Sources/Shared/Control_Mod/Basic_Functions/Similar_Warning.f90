!==============================================================================!
  subroutine Control_Mod_Similar_Warning(keyword, item, key_type)
!------------------------------------------------------------------------------!
!>  Checks if item in the argument list is similar to a keyword.  This was
!>  introduced to warn users that they might have misspelled a keyword in
!>  the control file. (It did not prove to be very useful though.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control   !! parent class
  character(len=*)           :: keyword   !! keyword being compared
  character(len=*)           :: item      !! item being compared
  character(len=*), optional :: key_type  !! type of the keyword
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
  logical :: found
!==============================================================================!

  ! Make monitoring points and porosity zones exceptions
  if(len_trim(keyword) .eq. 20) then
    if(keyword(1:16) .eq. 'MONITORING_POINT') return
  end if
  if(len_trim(keyword) .eq. 17) then
    if(keyword(1:13) .eq. 'POROUS_REGION') return
  end if

  !----------------------------------------------!
  !   If this item has already been recognized   !
  !     as "similar", do nothing, just return    !
  !----------------------------------------------!
  do n = 1, n_similar
    if( similar(n) .eq. item ) return
  end do

  !---------------------------------------------!
  !   This item is not in the list of similar   !
  !        ones, print a warning message        !
  !---------------------------------------------!
  if( Math % Approx_String(keyword, item, 1) ) then

    ! Print a warning message
    if(First_Proc()) then
      print *,  '#============================================='//           &
                '============================================='
      if(.not. present(key_type)) then
        print '(4a)', ' # NOTE! Could not find the keyword: ', keyword,      &
                      ', but found similar: ', item
      else
        print '(6a)', ' # NOTE! Could not find the ',key_type,': ',keyword,  &
                      ', but found similar: ', item
      end if
      print *, '#          Are you sure it is not a typing error'//          &
                   ' in the control file?'
      print *,  '#---------------------------------------------'//           &
                '---------------------------------------------'
    end if

    ! Store similar item in the list of similar items
    found = .false.
    do n = 1, n_similar
      if( similar(n) .eq. item ) found = .true.
    end do
    if(.not. found) then
      n_similar = n_similar + 1
      similar(n_similar) = item
    end if

  end if

  end subroutine
