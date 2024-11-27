!==============================================================================!
  subroutine Bulk_Fill(Info, Flow)
!------------------------------------------------------------------------------!
!>  Fills the previously created information box with actual bulk flow data
!>  from the simulation. It writes specific data like maximum Courant and
!>  Peclet numbers, volume flow rates in x, y, and z directions, and pressure
!>  drops in these directions as well. The data is formatted to fit into the
!>  structure created by Bulk_Start
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)             :: Info  !! parent, singleton object Info
  type(Field_Type), intent(in) :: Flow  !! flow field for which
                                        !! the info is fille
!==============================================================================!

  if(First_Proc()) then

    ! Courant and Peclet numbers
    write(Info % bulk % line(1)( 27: 49),    '(a23)') 'Maximum Courant number:'
    write(Info % bulk % line(1)( 51: 60), '(1pe9.3)') Flow % cfl_max

    write(Info % bulk % line(1)( 69: 90),    '(a22)') 'Maximum Peclet number:'
    write(Info % bulk % line(1)( 92:100), '(1pe9.3)') Flow % pe_max

    write(Info % bulk % line(2)( 27: 34),      '(a8)') 'Bulk U :'
    write(Info % bulk % line(2)( 36: 45), '(1pe10.3)') Flow % bulk % u
    write(Info % bulk % line(2)( 55: 62),      '(a8)') 'Bulk V :'
    write(Info % bulk % line(2)( 64: 75), '(1pe10.2)') Flow % bulk % v
    write(Info % bulk % line(2)( 83: 90),      '(a8)') 'Bulk W :'
    write(Info % bulk % line(2)( 92:101), '(1pe10.2)') Flow % bulk % w

    write(Info % bulk % line(3)( 27: 34),      '(a8)') 'Pdrop x:'
    write(Info % bulk % line(3)( 36: 45), '(1pe10.3)') Flow % bulk % p_drop_x
    write(Info % bulk % line(3)( 55: 62),      '(a8)') 'Pdrop y:'
    write(Info % bulk % line(3)( 64: 75), '(1pe10.2)') Flow % bulk % p_drop_y
    write(Info % bulk % line(3)( 83: 90),      '(a8)') 'Pdrop z:'
    write(Info % bulk % line(3)( 92:101), '(1pe10.2)') Flow % bulk % p_drop_z
  end if

  end subroutine
