#include "../../Shared/Assert.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Time_Mod
!------------------------------------------------------------------------------!
!>  This module manages the data and procedures related to time stepping in
!>  simulations, providing a centralized control mechanism for managing the
!>  progression of simulation time.  Module entails the definition of the
!>  Time_Type, global singleton object Time and a lot of time-step-specific
!>  data and methods.  The Time object is primarily used in the main function
!>  of the Process sub-program to manage time integration within the main
!>  solution loop, but due to its global accessibility, it can be utilized in
!>  other functions where decisions depend on the current iteration state.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Time type   !
  !---------------!
  !> A type that holds time-stepping-specific data and methods and is
  !> the basis for definition of the global singleton object Time
  type Time_Type

    real,    private :: physical_time      !! physical time of simulation [s]
    integer, private :: current_time_step  !! current dt in this simulation
    integer, private :: first_time_step    !! first dt in this simulation
    integer, private :: last_time_step     !! last dt in this simulation

    contains
      procedure :: Curr_Dt
      procedure :: First_Dt
      procedure :: Get_Time
      procedure :: Last_Dt
      procedure :: Needs_More_Steps
      procedure :: Increase_Time
      procedure :: Set_Curr_Dt
      procedure :: Set_First_Dt
      procedure :: Set_Last_Dt
      procedure :: Set_Time

  end type

  !----------------------!
  !   Singleton object   !
  !----------------------!
  type(Time_Type) :: Time  !! singleton type object of type Time_Type defined
                           !! globally for synchronization of its data members
                           !! and easier access to its member functions
  contains
#   include "Time_Mod/Curr_Dt.f90"
#   include "Time_Mod/First_Dt.f90"
#   include "Time_Mod/Get_Time.f90"
#   include "Time_Mod/Last_Dt.f90"
#   include "Time_Mod/Needs_More_Steps.f90"
#   include "Time_Mod/Increase_Time.f90"
#   include "Time_Mod/Set_Curr_Dt.f90"
#   include "Time_Mod/Set_First_Dt.f90"
#   include "Time_Mod/Set_Last_Dt.f90"
#   include "Time_Mod/Set_Time.f90"

  end module
