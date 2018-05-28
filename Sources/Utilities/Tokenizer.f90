!==============================================================================!
  module Tokenizer_Mod
!------------------------------------------------------------------------------!
!   This module is for tokenizing strings.                                     !
!   Useful when reading command files, user input or other input files.        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------------------!
  !   Boundary_Conditions type   !
  !------------------------------!
  type Tokenizer_Type
    character(len=300) :: whole              ! whole string
    character(len=300) :: tokens(300)        ! tokens             
    integer            :: n_tokens           ! number of tokens 
    integer            :: s(300), e(300)     ! tokens starts and ends
  end type

  type(Tokenizer_Type) :: line 
  integer, parameter   :: CMN_FILE = 11
  integer              :: cmn_line_count     ! command line count

  contains

!==============================================================================!
  subroutine Tokenizer_Mod_Read_Line(un, reached_end) 
!------------------------------------------------------------------------------!
!  Reads a line from a file (unit 9) and discards if it is comment.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: un  ! unit
  logical, optional :: reached_end
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  !-----------------------------------!
  !  Read the whole line into whole  !
  !-----------------------------------!
1 read(un,'(A300)', end=2) line % whole

  ! Shift the whole line to the left (remove leading spaces)
  line % whole = adjustl(line % whole)

  !------------------------------------------!
  !  If you are reading from command file    !
  !  (T-FlowS.cmn), increase the line count  !
  !------------------------------------------!
  if( un .eq. CMN_FILE ) then
    cmn_line_count = cmn_line_count + 1
  end if

  !--------------------!
  !  Skip empty lines  !
  !--------------------!
  if( line % whole  .eq.  '' ) goto 1 ! see: man ascii

  !----------------------!
  !  Skip comment lines  !
  !----------------------!
  if( trim(line % whole(1:1)) .eq. '!' .or.               &
      trim(line % whole(1:1)) .eq. '#' .or.               &
      trim(line % whole(1:1)) .eq. '%' ) goto 1

  !--------------------------------------!
  !  Parse tokens. This is somehow cool  !
  !--------------------------------------!
  line % n_tokens = 0
  if(line % whole(1:1) >= '!') then
    line % n_tokens = 1
    line % s(1)=1
  end if
  do i=1,298
    if( line % whole(i:  i  ) <  '!' .and.  &
        line % whole(i+1:i+1) >= '!') then
      line % n_tokens = line % n_tokens + 1
      line % s(line % n_tokens) = i+1
    end if
    if( line % whole(i  :i  ) >= '!' .and.  &
        line % whole(i+1:i+1) <  '!') then
      line % e(line % n_tokens) = i
    end if
  end do

  ! Chop them up
  do i = 1, line % n_tokens
    line % tokens(i) = line % whole(line % s(i) : line % e(i))
  end do

  return

  ! Error trap, if here, you reached the end of file
2 if(present(reached_end)) then
    reached_end = .true.
  end if

  end subroutine

  end module
