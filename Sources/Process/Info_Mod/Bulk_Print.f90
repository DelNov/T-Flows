!==============================================================================!
  subroutine Bulk_Print(Info, Flow, dom, n_dom)
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)             :: Info
  type(Field_Type), intent(in) :: Flow
  integer,          intent(in) :: dom, n_dom
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
