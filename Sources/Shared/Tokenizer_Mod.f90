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
  integer              :: cmn_line_count     ! command line count

  contains

  include 'Tokenizer_Mod/Read_Line.f90'

  end module
