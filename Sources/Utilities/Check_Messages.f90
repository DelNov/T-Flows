
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

  call Message % Print_Plain_Text(70,                                     &
                               "This is a message shown in the     "  //  &
                               "plain text mode.  Here, I will     "  //  &
                               "also show how new line is defined. "  //  &
                               "After this sentence, I will put    "  //  &
                               "a special sign for new line. \n    "  //  &
                               "I hope it did break the previous   "  //  &
                               "line.  If it didn't, now I will    "  //  &
                               "use three in a row. \n \n \n       "  //  &
                               "Now it should be obviuos.          "      &
                               )

  call Message % Print_Framed_Text(69, "",                                &
                               "This is a message without header,  "  //  &
                               "but with a frame around it.  As    "  //  &
                               "of now, the frame is simple, yet   "  //  &
                               "in accord with the rest of T-Flows "  //  &
                               "headers.  Thick line in the lead,  "  //  &
                               "thin line trailing.  If the header "  //  &
                               "is empy string, it will be skipped."      &
                               )

  call Message % Print_Framed_Text(59,                                    &
                               "THIS IS A MESSAGE WITH HEADER!",          &
                               "To get his look, one has to call   "  //  &
                               "the same function as before, but   "  //  &
                               "with non-empty argument for header "  //  &
                               "so that the function knows it has  "  //  &
                               "to be printed.  It also inserts a  "  //  &
                               "neat dash line in between the      "  //  &
                               "header and the rest of the text.   "  //  &
                               "This is one very long      message."      &
                               )

  call Message % Print_Warning(40,                                        &
                               "An example of a warning message.   "  //  &
                               "The philopshy is to call its       "  //  &
                               "function  without any parameters   "  //  &
                               "for header, but it insertst it on  "  //  &
                               "its own.  \n \n  I hope the new    "  //  &
                               "lines are still working and that   "  //  &
                               "you can see them above. With this  "  //  &
                               "function one can, however, use the "  //  &
                               "optional arguments in_file and     "  //  &
                               "at_line for easier debugging. But, "  //  &
                               "macros __FILE___ and __LINE__ do   "  //  &
                               "not seem to work for included files",     &
                               in_file="Main_Con.f90", at_line=53)

  call Message % Print_Error  (60,                                        &
                               "This is an error long message. It  "  //  &
                               "is invoked in the same way as her  "  //  &
                               "warning sister, but an important   "  //  &
                               "difference between them is that    "  //  &
                               "the warning sister does not stop   "  //  &
                               "the execution of the program, but  "  //  &
                               "this one does. \n \n For a descent "  //  &
                               "exit, one should uses the MPI      "  //  &
                               "calls to end parallel execution,   "  //  &
                               "followed by a stop statement.",           &
                               in_file="Main_Con.f90", at_line=68)

  end program
