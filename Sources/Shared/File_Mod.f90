#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module File_Mod
!------------------------------------------------------------------------------!
!>  This module contains variables, subroutines and function for managing file
!>  operations, including reading, writing, deleting, and modifying files in
!>  both ASCII and binary formats.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use String_Mod
  use Message_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITEMS   =    2048  !! maximum items in work arrays
  integer, parameter :: BUFFER_SIZE = 1048576  !! buffer size
!==============================================================================!

  character(SL) :: problem_name(MD)  !! stores problem names for all domains

  !---------------!
  !   File type   !
  !---------------!
  !> File_Type encapsulates a set of procedures for various file
  !> operations and also buffers (i_buffer for integers, r_buffer
  !> for reals) used in buffered I/O.
  type File_Type

    ! Buffers for reading
    integer, private :: i_buffer(BUFFER_SIZE)  !! integer buffer
    real,    private :: r_buffer(BUFFER_SIZE)  !! real buffer

    contains
      procedure :: Append_For_Writing_Ascii
      procedure :: Buffered_Read_Int_Array
      procedure :: Buffered_Read_Real_Array
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
      procedure :: Single_Int_From_Keyboard
      procedure :: Single_Word_From_Keyboard

  end type

  !---------------------------!
  !   Singleton object File   !
  !---------------------------!
  type(File_Type) :: File  !! singleton File object allowing easy access to
                           !! module's functionalities across the application

  integer(SP) :: int4_array(MAX_ITEMS)   !! working integer(4) array
  integer(DP) :: int8_array(MAX_ITEMS)   !! working integer(8) array
  real(SP)    :: real4_array(MAX_ITEMS)  !! working real(4) array
  real(DP)    :: real8_array(MAX_ITEMS)  !! working real(8) array

  contains

#   include "File_Mod/Append_For_Writing_Ascii.f90"
#   include "File_Mod/Buffered_Read_Int_Array.f90"
#   include "File_Mod/Buffered_Read_Real_Array.f90"
#   include "File_Mod/Delete.f90"
#   include "File_Mod/Is_In_Binary_Format.f90"
#   include "File_Mod/Is_In_Unix_Format.f90"
#   include "File_Mod/Line_Length.f90"
#   include "File_Mod/Open_For_Reading_Ascii.f90"
#   include "File_Mod/Open_For_Reading_Binary.f90"
#   include "File_Mod/Open_For_Writing_Ascii.f90"
#   include "File_Mod/Open_For_Writing_Binary.f90"
#   include "File_Mod/Read_Binary_Int4_Array.f90"
#   include "File_Mod/Read_Binary_Int8_Array.f90"
#   include "File_Mod/Read_Binary_Real4_Array.f90"
#   include "File_Mod/Read_Binary_Real8_Array.f90"
#   include "File_Mod/Read_Line.f90"
#   include "File_Mod/Set_Name.f90"
#   include "File_Mod/Single_Int_From_Keyboard.f90"
#   include "File_Mod/Single_Word_From_Keyboard.f90"

  end module
