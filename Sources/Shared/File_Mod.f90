!==============================================================================!
  module File_Mod
!------------------------------------------------------------------------------!
!   Variables and subroutines for handling file access.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  character(len=80) :: problem_name(MD)

  !--------------------!
  !   Tokenizer type   !
  !--------------------!
  type Tokenizer_Type
    character(len=300) :: whole              ! whole string
    character(len=300) :: tokens(300)        ! tokens
    integer            :: n_tokens           ! number of tokens
    integer            :: s(300), e(300)     ! tokens starts and ends
  end type

  type(Tokenizer_Type) :: line

  contains

  include 'File_Mod/Append_File_For_Writing.f90'
  include 'File_Mod/Set_Name.f90'
  include 'File_Mod/Open_File_For_Reading.f90'
  include 'File_Mod/Open_File_For_Reading_Binary.f90'
  include 'File_Mod/Open_File_For_Writing.f90'
  include 'File_Mod/Open_File_For_Writing_Binary.f90'
  include 'File_Mod/Read_Line.f90'

  end module
