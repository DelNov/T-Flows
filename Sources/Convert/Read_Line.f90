!==============================================================================!
  subroutine Read_Line(un, string, tn, ts, te) 
!------------------------------------------------------------------------------!
!   Reads a line from a file (unit un) and parses it.                          ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Parameters]---------------------------------!
  integer   :: un
  character :: string*300
  integer   :: tn, ts(300), te(300)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  read(un,'(A300)') string

  ! Parse tokens. This is somehow cool 
  tn = 0
  if(string(1:1) >= '!') then
    tn = 1
    ts(1)=1
  end if
  do i=1,298
    if( string(i:i)  < '!' .and. string(i+1:i+1) >= '!') then
      tn = tn + 1
      ts(tn) = i+1
    end if
    if( string(i:i) >= '!' .and. string(i+1:i+1)  < '!') then
      te(tn) = i
    end if
  end do

  end subroutine Read_Line
