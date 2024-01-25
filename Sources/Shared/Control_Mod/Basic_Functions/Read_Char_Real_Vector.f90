!==============================================================================!
  subroutine Read_Char_Real_Vector(Control, keyword, n, defc, defv,  &
                                   valc, valv, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read character "valc" and vector "valv" behind a !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   value specified in argument "def" is used.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)              :: Control
  character(len=*),    intent(in)  :: keyword
  integer,             intent(in)  :: n        ! size of array (typically small)
  character(len=*),    intent(in)  :: defc     ! default value for character
  real,                intent(in)  :: defv(n)  ! default value for vector
  character(SL),       intent(out) :: valc     ! spefified value, if found
  real,                intent(out) :: valv(n)  ! spefified value, if found
  logical,   optional, intent(in)  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: reached_end
  integer :: i
!==============================================================================!

  rewind(Control % file_unit)

  ! Set default value
  valc = defc
  valv = defv

  !--------------------------------------------------------------!
  !   Browse through command file to see if user specificed it   !
  !--------------------------------------------------------------!
  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) goto 1

    ! Found the correct keyword
    if(Line % tokens(1) .eq. trim(keyword)) then
      read(Line % tokens(2), *) valc
      do i = 1, n
        read(Line % tokens(i+2), *, end=2, err=2) valv(i)
      end do
2     continue
      return
    end if

  end do

  !--------------------------------------------!
  !   Keyword was not found; issue a warning   !
  !--------------------------------------------!
1 if(present(verbose)) then
    if(verbose .and. First_Proc()) then
      print '(4a)', ' # NOTE! Could not find the keyword: ',  &
                      trim(keyword), '. Using the default: ', trim(defc)
    end if
  end if

  end subroutine
