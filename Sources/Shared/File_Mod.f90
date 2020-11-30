!==============================================================================!
  module File_Mod
!------------------------------------------------------------------------------!
!   Variables and subroutines for handling file access.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_TOKENS = DL / 2
!==============================================================================!

  character(SL) :: problem_name(MD)

  !--------------------!
  !   Tokenizer type   !
  !--------------------!
  type Tokenizer_Type
    character(DL) :: whole               ! whole string
    character(SL) :: tokens(MAX_TOKENS)  ! tokens
    integer       :: n_tokens            ! number of tokens
    integer       :: s(MAX_TOKENS),  &   ! tokens starts ...
                     e(MAX_TOKENS)       ! ... and ends
    character(1)  :: first, last         ! first and last characters in whole
  end type

  type(Tokenizer_Type) :: line

  contains

  include 'File_Mod/Append_File_For_Writing.f90'
  include 'File_Mod/Delete.f90'
  include 'File_Mod/Is_File_In_Unix_Format.f90'
  include 'File_Mod/Set_Name.f90'
  include 'File_Mod/Open_File_For_Reading.f90'
  include 'File_Mod/Open_File_For_Reading_Binary.f90'
  include 'File_Mod/Open_File_For_Writing.f90'
  include 'File_Mod/Open_File_For_Writing_Binary.f90'
  include 'File_Mod/Read_Line.f90'

  end module
