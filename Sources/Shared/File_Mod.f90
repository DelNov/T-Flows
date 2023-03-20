#include "../Shared/Unused.h90"

!==============================================================================!
  module File_Mod
!------------------------------------------------------------------------------!
!   Variables and subroutines for handling file access.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use String_Mod
  use Message_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITEMS = 2048
!==============================================================================!

  character(SL) :: problem_name(MD)

  !---------------!
  !   File type   !
  !---------------!
  type File_Type

    contains
      procedure :: Append_For_Writing_Ascii
      procedure :: Delete
      procedure :: Is_In_Binary_Format
      procedure :: Is_In_Unix_Format
      procedure :: Line_Length
      procedure :: Open_For_Reading_Ascii
      procedure :: Open_For_Reading_Binary
      procedure :: Open_For_Writing_Ascii
      procedure :: Open_For_Writing_Binary
      procedure :: Read_Binary_Int4_Array
      procedure :: Read_Binary_Int8_Array
      procedure :: Read_Binary_Real4_Array
      procedure :: Read_Binary_Real8_Array
      procedure :: Read_Line
      procedure :: Set_Name

  end type

  type(Tokenizer_Type) :: Line
  type(File_Type)      :: File

  integer(SP) :: int4_array(MAX_ITEMS)
  integer(DP) :: int8_array(MAX_ITEMS)
  real(SP)    :: real4_array(MAX_ITEMS)
  real(DP)    :: real8_array(MAX_ITEMS)

  contains

#   include "File_Mod/Append_For_Writing_Ascii.f90"
#   include "File_Mod/Delete.f90"
#   include "File_Mod/Is_In_Binary_Format.f90"
#   include "File_Mod/Is_In_Unix_Format.f90"
#   include "File_Mod/Line_Length.f90"
#   include "File_Mod/Set_Name.f90"
#   include "File_Mod/Open_For_Reading_Ascii.f90"
#   include "File_Mod/Open_For_Reading_Binary.f90"
#   include "File_Mod/Open_For_Writing_Ascii.f90"
#   include "File_Mod/Open_For_Writing_Binary.f90"
#   include "File_Mod/Read_Binary_Int4_Array.f90"
#   include "File_Mod/Read_Binary_Int8_Array.f90"
#   include "File_Mod/Read_Binary_Real4_Array.f90"
#   include "File_Mod/Read_Binary_Real8_Array.f90"
#   include "File_Mod/Read_Line.f90"

  end module
