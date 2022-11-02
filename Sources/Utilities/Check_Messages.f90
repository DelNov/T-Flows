
!==============================================================================!
  program Check_Messages
!------------------------------------------------------------------------------!
!   Compile as:                                                                !
!                                                                              !
!   gfortran -c -cpp ../Shared/Const_Mod.f90                                   !
!   gfortran -c -cpp ../Shared/Tokenizer_Mod.f90                               !
!   gfortran -c -cpp ../Shared/Comm_Mod_Seq.f90                                !
!   gfortran -c -cpp ../Shared/Message_Mod.f90                                 !
!   gfortran -o check_msg -cpp *.o Check_Messages.f90                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Message_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call Message % Framed(42,                                         &
                        "This message is 42 characters wide ",      &
                        "0 1 2 3 4 5 6 7 8 9                "   //  &
                        "0 1 2 3 4 5 6 7 8 9                "   //  &
                        "0 1 2 3 4 5 6 7 8 9                ")

  call Message % Framed(22,                                        &
                       "This message is only 22 characters "   //  &
                       "wide, but the header is wider and  "   //  &
                       "the message width is adapted for it",      &
                       "0 1 2 3 4 5 6 7 8 9                "   //  &
                       "0 1 2 3 4 5 6 7 8 9                "   //  &
                       "0 1 2 3 4 5 6 7 8 9                ")

  call Message % Frameless(70,                                       &
                           "This is a message shown in the     "  //  &
                           "frameless mode.   Here, I will     "  //  &
                           "also show how new line is defined. "  //  &
                           "After this sentence, I will put    "  //  &
                           "a special sign for new line. \n    "  //  &
                           "I hope it did break the previous   "  //  &
                           "line.  If it didn't, now I will    "  //  &
                           "use three in a row. \n \n \n       "  //  &
                           "Now it should be obviuos.          "      &
                           )

  call Message % Framed(69, "",                                    &
                        "This is a message without header,  "  //  &
                        "but with a frame around it.  As    "  //  &
                        "of now, the frame is simple, yet   "  //  &
                        "in accord with the rest of T-Flows "  //  &
                        "headers.  Thick line in the lead,  "  //  &
                        "thin line trailing.  If the header "  //  &
                        "is empy string, it will be skipped."      &
                        )

  call Message % Framed(59,                                        &
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

  call Message % Warning(40,                                        &
                         "An example of a warning message.   "  //  &
                         "The philopshy is to call its       "  //  &
                         "function  without any parameters   "  //  &
                         "for header, but it insertst it on  "  //  &
                         "its own.  \n \n  I hope the new    "  //  &
                         "lines are still working and that   "  //  &
                         "you can see them above. With this  "  //  &
                         "function one can, however, use the "  //  &
                         "optional arguments in_file and     "  //  &
                         "at_line for easier debugging.      ",     &
                         file=__FILE__, line=__LINE__)

  call Message % Error(60,                                        &
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
                       file=__FILE__, line=__LINE__)

  end program
