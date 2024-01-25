!==============================================================================!
  subroutine Mass_Transfer(Control, mt_model, c_lee, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer              :: mt_model  ! 0: no mass transfer, 1 for temperature gradient, 2 for Lee
  real                 :: c_lee(2)
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: valc
  real          :: defv(2), valv(2)
!==============================================================================!

  ! set initial values
  mt_model = 0
  data defv / -1.0, -1.0 /

  !call Control % Read_Char_Item('MASS_TRANSFER', 'no', val, verbose)
  call Control % Read_Char_Real_Vector('MASS_TRANSFER', 2, 'no', defv,  &
                                      valc, valv, verbose)
  call String % To_Upper_Case(valc)

  if( valc .eq. 'YES' ) then
    mt_model = 1  ! temperature gradient model
    if(First_Proc())  &
      print '(a)', ' # Mass_Transfer: m_dot is based on temperature gradient'

  else if( valc .eq. 'LEE' ) then
    mt_model = 2  ! Lee
    if (valv(1)<0.0 .or. valv(2)<0.0) then
      call Message % Error(60,                                               &
          'Lee model is used for MASS_TRANSFER but its coefficients are '//  &
          'not given. Two coefficients are necessary:'//           &
          'e.g. MASS_TRANSFER LEE 100 1000 '//                     &
          'This error is critical.  Exiting.',                               &
          file=__FILE__, line=__LINE__, one_proc=.true.)
    endif
    c_lee(1) = valv(1)
    c_lee(2) = valv(2)
    if(First_Proc()) then
      print '(a,PE12.4)', ' # Mass_Transfer: C_Lee for condensation= ', c_lee(1)
      print '(a,PE12.4)', ' # Mass_Transfer: C_Lee for vaporization= ', c_lee(2)
    endif

  else if( valc .eq. 'NO' ) then
    mt_model = 0

  else
    call Message % Error(60,                                              &
                         'Unknown state for phase change: '//trim(valc)// &
                         '. \n This error is critical.  Exiting.',        &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine

