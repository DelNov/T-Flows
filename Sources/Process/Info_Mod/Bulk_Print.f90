!==============================================================================!
  subroutine Bulk_Print(Info, Flow, dom, n_dom)
!------------------------------------------------------------------------------!
!>  Prints the filled information box on the screen. It ensures that the
!>  formatted data, including the bulk flow parameters and numerical values,
!>  are displayed correctly. This subroutine may be called many times during
!>  the simulation to update and display the latest bulk flow information.
!>  The printing is conditional on the domain of the simulation (indicated by
!>  dom and n_dom), allowing selective information display for each domain
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)             :: Info   !! parent, singleton object Info
  type(Field_Type), intent(in) :: Flow   !! flow field for info is printed
  integer,          intent(in) :: dom    !! domain for which info is printed
  integer,          intent(in) :: n_dom  !! number of domains in the simulation
!==============================================================================!

  call Info % Bulk_Fill(Flow)

  if(First_Proc()) then

    ! String is L_LINE+2 long
    if(dom .eq. 1) then
      print '(a129)', trim(Info % bulk % line_lead)
    else
      print '(a108)', trim(Info % bulk % line_foll)
    end if
    print '(a108)', trim(Info % bulk % line(1))
    print '(a108)', trim(Info % bulk % line_sep)
    print '(a108)', trim(Info % bulk % line(2))
    print '(a108)', trim(Info % bulk % line(3))
    if(dom .eq. n_dom) then
      print '(a108)', trim(Info % bulk % line_trail)
      print '(a)',    ''
    end if

  end if

  end subroutine
