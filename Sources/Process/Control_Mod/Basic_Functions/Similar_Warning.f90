!==============================================================================!
  subroutine Control_Mod_Similar_Warning(keyword, item, verbose, key_type)
!------------------------------------------------------------------------------!
!   Checks if item in the argument list is similar to a keyword.               !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Comm_Mod, only: this_proc
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)           :: keyword
  character(len=*)           :: item
  logical,          optional :: verbose
  character(len=*), optional :: key_type
!---------------------------------[Interface]----------------------------------!
  include '../Shared/Approx_String.int'
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
  logical :: found
!==============================================================================!

  ! Make monitoring points an exception
  if(len_trim(keyword) .eq. 20) then
    if(keyword(1:16) .eq. 'MONITORING_POINT') return
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
  if( Approx_String(keyword, item, 1) ) then

    ! Print a warning message
    if(this_proc < 2) then
      print '(a)',  ' #============================================='//      &
                    '============================================='
      if(.not. present(key_type)) then
        print '(4a)', ' # NOTE! Could not find the keyword: ', keyword,      &
                      ', but found similar: ', item
      else
        print '(6a)', ' # NOTE! Could not find the ',key_type,': ',keyword,  &
                      ', but found similar: ', item
      end if
      print '(a)', ' #          Are you sure it is not a typing error'//     &
                   ' in the control file?'
      print '(a)',  ' #---------------------------------------------'//      &
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
