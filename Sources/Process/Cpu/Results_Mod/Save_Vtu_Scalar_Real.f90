!==============================================================================!
  subroutine Save_Vtu_Scalar_Real(Results,                &
                                  var_name, plot_inside,  &
                                  val, fs, fp,            &
                                  data_offset, sweep)
!------------------------------------------------------------------------------!
!>  The Save_Vtu_Scalar_Real subroutine is designed to write a single real
!>  scalar variable defined over cells into a .vtu file as a single record.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparatory work: Connects a real cell buffer for data handling and      !
!     sets the plotting precision.                                             !
!   * Data handling: Based on the 'sweep' stage, it either writes the header   !
!     data for the variable or the actual scalar values to the file.           !
!   * Data writing: Scalar values are written into the file in an appended     !
!     format, taking into account whether the data is for internal or boundary !
!     cells.                                                                   !
!   * Data offset management: Updates the data offset after writing, ensuring  !
!     accurate placement of subsequent data in the file.                       !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * Buffer connection: A buffer is connected to handle the scalar values.    !
!   * Header writing: In the first sweep, the subroutine writes the XML header !
!     for the scalar variable, including its name and data type.               !
!   * Data writing: In the second sweep, it writes the actual scalar values.   !
!     This includes processing values for either internal or boundary cells    !
!     based on the 'plot_inside' flag.                                         !
!   * Closing operations: Finally, it updates the data offset and disconnects  !
!     the buffer, completing the subroutine's operations.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results      !! parent class
  character(len=*)    :: var_name     !! name of the variable
  logical             :: plot_inside  !! true to plots inside,
                                      !! false to plot on the boundary
  real                :: val(:)       !! variable's values
  integer             :: fs, fp       !! file unit sequential and parallel
  integer             :: data_offset  !! offset in the .vtu file
  integer             :: sweep        !! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)               :: data_size
  integer                   :: i, c1, c2, c_f, c_l
  character(SL)             :: str1
  real, pointer, contiguous :: buffer(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  call Work % Connect_Real_Cell(buffer)

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val, 1)
  c_l = ubound(val, 1)

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
      write(fs) IN_4                                  //  &
                '<PDataArray type='//floatp           //  &
                ' Name="' // trim(var_name) // '"/>'  // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//floatp         //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then

    call Profiler % Start('Save_Vtu_Results (scalar real)')

    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP, SP)
      write(fp) data_size
      i = 0
      do c1 = c_f, c_l
        i = i + 1
        buffer(i) = val(c1)
      end do
      write(fp) buffer(1:i)
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP, SP)
      end do
      write(fp) data_size
      i = 0
      do c2 = c_f, c_l
        i = i + 1
        buffer(i) = val(c2)
      end do
      write(fp) buffer(1:i)
    end if

    call Profiler % Stop('Save_Vtu_Results (scalar real)')

  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * RP
    else
      do c2 = c_f, c_l
        data_offset = data_offset + RP
      end do
    end if
    data_offset = data_offset + SP
  end if

  call Work % Disconnect_Real_Cell(buffer)

  end subroutine
