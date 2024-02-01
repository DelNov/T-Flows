!==============================================================================!
  subroutine Save_Vtu_Vector_Real(Results,                      &
                                  var_name, plot_inside,        &
                                  val_1, val_2, val_3, fs, fp,  &
                                  data_offset, sweep)
!------------------------------------------------------------------------------!
!>  Save_Vtu_Vector_Real is engineered to write a real vector variable defined
!>  over cells into a .vtu file as a single record.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparatory work: It sets up precision for plotting and allocates a      !
!     buffer for data handling.                                                !
!   * Data processing: Based on the sweep stage, it writes either the XML      !
!     header for the vector variable or the actual vector data to the file.    !
!   * Data writing: Vector values are written in an appended format,           !
!     considering whether the data is for inside or boundary cells.            !
!   * Data offset update: Updates the data offset post writing, ensuring       !
!     correct placement of subsequent data in the file.                        !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * Buffer allocation: Allocates a buffer to store the vector data.          !
!   * Header writing: In the first sweep, the subroutine writes the XML header !
!     for the vector variable, including its name and data format.             !
!   * Data writing: In the second sweep, it writes the actual vector values.   !
!     This step accounts for whether the data is for internal or boundary cells!
!     based on the 'plot_inside' flag.                                         !
!   * Closing operations: Finally, updates the data offset and completes the   !
!     subroutine's tasks.                                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results      !! parent class
  character(len=*)    :: var_name     !! name of the variable
  logical             :: plot_inside  !! true to plots inside,
                                      !! false to plot on the boundary
  real                :: val_1(:)     !! vector component
  real                :: val_2(:)     !! vector component
  real                :: val_3(:)     !! vector component
  integer             :: fs, fp       !! file unit sequential and parallel
  integer             :: data_offset  !! offset in .vtu file
  integer             :: sweep        !!! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c1, c2, c_f, c_l
  character(SL) :: str1
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val_1, 1)
  c_l = ubound(val_1, 1)

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
      write(fs) IN_4                                //  &
                '<PDataArray type='//floatp         //  &
                ' Name="' // trim(var_name) // '"'  //  &
                ' NumberOfComponents="3"/>'         // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//floatp         //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' NumberOfComponents="3"'          //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP * 3, SP)
      write(fp) data_size
      do c1 = c_f, c_l
        write(fp) val_1(c1), val_2(c1), val_3(c1)
      end do
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP * 3, SP)
      end do
      write(fp) data_size
      do c2 = c_f, c_l
        write(fp) val_1(c2), val_2(c2), val_3(c2)
      end do
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * RP * 3
    else
      do c2 = c_f, c_l
        data_offset = data_offset + RP * 3
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine
