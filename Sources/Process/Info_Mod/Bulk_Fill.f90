!==============================================================================!
  subroutine Info_Mod_Bulk_Fill(Flow)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: Flow
!==============================================================================!

  if(First_Proc()) then

    ! Courant and Peclet numbers
    write(bulk_info % lines(1)( 27: 49),    '(a23)') 'Maximum Courant number:'
    write(bulk_info % lines(1)( 51: 60), '(1pe9.3)') Flow % cfl_max

    write(bulk_info % lines(1)( 69: 90),    '(a22)') 'Maximum Peclet number:'
    write(bulk_info % lines(1)( 92:100), '(1pe9.3)') Flow % pe_max

    write(bulk_info % lines(2)( 27: 34),      '(a8)') 'Flux x :'
    write(bulk_info % lines(2)( 36: 45), '(1pe10.3)') Flow % bulk % flux_x
    write(bulk_info % lines(2)( 55: 62),      '(a8)') 'Flux y :'
    write(bulk_info % lines(2)( 64: 75), '(1pe10.2)') Flow % bulk % flux_y
    write(bulk_info % lines(2)( 83: 90),      '(a8)') 'Flux z :'
    write(bulk_info % lines(2)( 92:101), '(1pe10.2)') Flow % bulk % flux_z

    write(bulk_info % lines(3)( 27: 34),      '(a8)') 'Pdrop x:'
    write(bulk_info % lines(3)( 36: 45), '(1pe10.3)') Flow % bulk % p_drop_x
    write(bulk_info % lines(3)( 55: 62),      '(a8)') 'Pdrop y:'
    write(bulk_info % lines(3)( 64: 75), '(1pe10.2)') Flow % bulk % p_drop_y
    write(bulk_info % lines(3)( 83: 90),      '(a8)') 'Pdrop z:'
    write(bulk_info % lines(3)( 92:101), '(1pe10.2)') Flow % bulk % p_drop_z
  end if

  end subroutine
