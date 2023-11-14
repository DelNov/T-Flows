#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Message_Mod
!------------------------------------------------------------------------------!
!>  This module is designed to provide a suite of procedures for handling and
!>  displaying messages in a structured and formatted manner. It is focused on
!>  printing warning and error messages, as well as creating framed messages
!>  with different line styles.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Tokenizer_Mod
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Message type   !
  !------------------!
  !> Message_Type includes several procedures for different
  !> aspects of message handling and formatting.
  type Message_Type

    contains
      procedure, private :: Dashed_Line
      procedure          :: Error
      procedure          :: Framed
      procedure, private :: Frameless
      procedure          :: Warning
      procedure, private :: Thick_Line
      procedure, private :: Thin_Line

  end type

  type(Message_Type) :: Message  !! A global object of the type Message_Type

  contains

#   include "Message_Mod/Dashed_Line.f90"
#   include "Message_Mod/Error.f90"
#   include "Message_Mod/Framed.f90"
#   include "Message_Mod/Frameless.f90"
#   include "Message_Mod/Warning.f90"
#   include "Message_Mod/Thick_Line.f90"
#   include "Message_Mod/Thin_Line.f90"

  end module
