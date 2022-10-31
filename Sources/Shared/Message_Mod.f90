!==============================================================================!
  module Message_Mod
!------------------------------------------------------------------------------!
!   Procedures for printing warning and error messages.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Tokenizer_Mod
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Message type   !
  !------------------!
  type Message_Type

    contains
      procedure, private :: Dashed_Line
      procedure          :: Print_Error
      procedure          :: Print_Framed_Text
      procedure          :: Print_Plain_Text
      procedure          :: Print_Warning
      procedure, private :: Thick_Line
      procedure, private :: Thin_Line

  end type

  type(Message_Type) :: Message

  contains

  include 'Message_Mod/Dashed_Line.f90'
  include 'Message_Mod/Print_Error.f90'
  include 'Message_Mod/Print_Framed_Text.f90'
  include 'Message_Mod/Print_Plain_Text.f90'
  include 'Message_Mod/Print_Warning.f90'
  include 'Message_Mod/Thick_Line.f90'
  include 'Message_Mod/Thin_Line.f90'

  end module
