!==============================================================================!
  subroutine User_Mod_Get_User_Field_For_Saving(Grid,                     &
                                                rank,                     &
                                                c_f, c_l, field_to_save,  &
                                                field_name)
!------------------------------------------------------------------------------!
!   This function is called while save the fields in .vtu file format          !
!                                                                              !
!   Note: Be aware that this function is called from CPU only, at the moment   !
!         when results are downloaded to host (CPU).  Thus, make sure that     !
!         you send the CPU copy of the variable.                    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
  integer         :: rank                    ! user variable rank (1 - 60)
  integer         :: c_f, c_l
  real            :: field_to_save(c_f:c_l)
  character(SL)   :: field_name
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  !------------------------------------------------------------!
  !                                                            !
  !   Here you can choose up to 60 user variables for saving   !
  !                                                            !
  !------------------------------------------------------------!
  if(rank .eq. 1) then

    ! Keep in mind that user fields must be defined
    ! in the range -Grid % n_bnd_cells:Grid % n_cells
    !@ Assert(lbound(my_useful_field,1) .eq. -Grid % n_bnd_cells)
    !@ Assert(ubound(my_useful_field,1) .eq.  Grid % n_cells)

    ! Copy user field to Results_Mod's memory space
    !@ do c = c_f, c_l
    !@   field_to_save(c) = my_useful_field(c)
    !@ end do

    ! Set the user field's name - this will appear in .vtu file
    ! If the line setting "field_name" is ommited, the field won't
    ! be saved.  It is commented by default, but if you want your
    ! variable to be saved, uncomment it.
    ! field_name = "My Useful Field [m^2/s^2]"

  end if

  end subroutine
