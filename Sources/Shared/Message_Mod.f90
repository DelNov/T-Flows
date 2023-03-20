#include "../Shared/Assert.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Message_Mod
!------------------------------------------------------------------------------!
!   Procedures for printing warning and error messages.                        !
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

  type(Message_Type) :: Message

  contains

#   include "Message_Mod/Dashed_Line.f90"
#   include "Message_Mod/Error.f90"
#   include "Message_Mod/Framed.f90"
#   include "Message_Mod/Frameless.f90"
#   include "Message_Mod/Warning.f90"
#   include "Message_Mod/Thick_Line.f90"
#   include "Message_Mod/Thin_Line.f90"

  end module
