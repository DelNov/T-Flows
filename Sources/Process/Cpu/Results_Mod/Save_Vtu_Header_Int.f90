!==============================================================================!
  subroutine Save_Vtu_Header_Int(Results, var_name, fs, fp, data_offset)
!------------------------------------------------------------------------------!
!>  The Save_Vtu_Header_Int subroutine is designed to write a header for
!>  and integer scalar variable into a .vtu and .pvtu files.  This is useful
!>  for variables which are not cell based, which should not rely on updates
!>  of data_offset as it is done in Save_Vtu_Scalar_Int subroutine.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results      !! parent class
  character(len=*)    :: var_name     !! name of the variabl
  integer             :: fs, fp       !! file unit sequential and parallel
  integer             :: data_offset  !! data offset in the .vtu file
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: str1
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  if(Parallel_Run() .and. First_Proc()) then
    write(fs) IN_4                                  //  &
              '<PDataArray type='//intp             //  &
              ' Name="' // trim(var_name) // '"/>'  // LF
  end if

  write(str1, '(i0.0)') data_offset
  if(data_offset .eq. 0) write(str1, '(i1)') data_offset
  write(fp) IN_4 // '<DataArray type='//intp           //  &
                    ' Name="' // trim(var_name) // '"' //  &
                    ' format="appended"'               //  &
                    ' offset="' // trim(str1) //'">'   // LF
  write(fp) IN_4 // '</DataArray>' // LF

  end subroutine
