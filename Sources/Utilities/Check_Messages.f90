
!==============================================================================!
  program Check_Messages
!------------------------------------------------------------------------------!
!   Compile as:                                                                !
!                                                                              !
!   gfortran -c ../Shared/Const_Mod.f90                                        !
!   gfortran -c ../Shared/Tokenizer_Mod.f90                                    !
!   gfortran -c ../Shared/Comm_Mod_Seq.f90                                     !
!   gfortran -c ../Shared/Message_Mod.f90                                      !
!   gfortran -o check_msg *.o Check_Messages.f90                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Message_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call Message % Print_Plain_Text(39,                                     &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "   //  &
                               "This is a box with plain text ... "       &
                               )

  call Message % Print_Framed_Text(69, "",                                &
                               "This is a message without header. "   //  &
                               "This is a message without header. "   //  &
                               "This is a message without header. "   //  &
                               "This is a message without header. "   //  &
                               "This is a message without header. "   //  &
                               "This is a message without header. "   //  &
                               "This is a message without header. "       &
                               )

  call Message % Print_Framed_Text(69,                                 &
                               "THIS IS A MESSAGE WITH HEADER!",       &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message. "   //  &
                               "This is one very long message."        &
                               )

  call Message % Print_Warning(40,                                       &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message. "   //  &
                               "This is an warning long message.",       &
                               in_file="Main_Con.f90", at_line=53)

  call Message % Print_Error  (60,                                     &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message. "   //  &
                               "This is an error long message.",       &
                               in_file="Main_Con.f90", at_line=53)

  end program
