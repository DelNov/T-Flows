!==============================================================================!
  subroutine Create_Porosity(Por, Grid)
!---------------------------------[Arguments]----------------------------------!
  class(Porosity_Type)    :: Por
  type(Grid_Type), target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: reg, n, c, fu
  logical       :: found
  character(SL) :: porous_region_rank
  character(SL) :: next_strings(MAX_STRING_ITEMS)
!==============================================================================!

  !---------------------------!
  !   Compiled from Convert   !  Note: T_FLOWS_PROGRAM 1 is Convert
  !---------------------------!
# if T_FLOWS_PROGRAM == 1
  print '(a)', " #===================================================="
  print '(a)', " # You are dealing with a problem with porous regions "
  print '(a)', " #----------------------------------------------------"
  print '(a)', " # Enter the number of porous regions: "
  print '(a)', " #----------------------------------------------------"
  Por % n_regions = File % Single_Int_From_Keyboard()

  !---------------------------!
  !   Compiled from Process   !  Note: T_FLOWS_PROGRAM 4 is Process
  !---------------------------!
# elif T_FLOWS_PROGRAM == 4
  ! Read number of Porous regions from control file
  call Control % Read_Int_Item('NUMBER_OF_POROUS_REGIONS', 0, &
                               Por % n_regions, .true.)
# endif

  if(Por % n_regions .eq. 0) return

  Por % pnt_grid => Grid

  !-----------------------------------------------------!
  !   Allocate memory (for both Convert and Processs)   !
  !-----------------------------------------------------!
  allocate(Por % region(Por % n_regions))

  ! Mahir's correction
  do reg = 1, Por % n_regions
    allocate(Por % region(reg) % cell_porous(Grid % n_cells))
    Por % region(reg) % cell_porous(:) = .false.
  end do

  !----------------------------------------!
  !                                        !
  !   Read porous regions characterisics   !
  !                                        !
  !----------------------------------------!
# if T_FLOWS_PROGRAM == 1
  do reg = 1, Por % n_regions
    print '(a,i2)', " # Enter the name of the STL file for porous region ", reg
    Por % region(reg) % stl_name = File % Single_Word_From_Keyboard()
  end do
# endif

# if T_FLOWS_PROGRAM == 4
  do reg = 1, Por % n_regions

    ! Set region's name
    write(porous_region_rank, '(a,i3.3)') 'POROUS_REGION_', reg

    ! Look for it
    call Control % Position_At_One_Key(porous_region_rank, found, .true.)

    !------------------------------------------!
    !   Found the section with porous region   !
    !------------------------------------------!
    if(found) then

      call Control % Read_Keyless_Strings_On(next_strings, n, .true.)

      call String % To_Upper_Case(next_strings(1))
      if(next_strings(1) .eq. 'STL_FILE') then
        Por % region(reg) % stl_name = next_strings(2)
        print '(a)', trim(Por % region(reg) % stl_name)
      else if(next_strings(1) .eq. 'FROM_GRID') then
        Por % region(reg) % stl_name = ''           ! empty, undefined
        ! Do some checks here
        found = any(Grid % por(:) .eq. reg)
        call Global % Lor_Log(found)  ! any processor found it?
        if(.not. found) then          ! not good for parallel
          call Message % Error(72,                                      &
                 'You specified keyword "FROM_GRID" in the porous ' //  &
                 'regions of the control file, but the grid does  ' //  &
                 'not have information on porous regions. \n \n   ' //  &
                 'This error is critical.  Exiting!',                   &
                 file=__FILE__, line=__LINE__)
        end if
      end if

    else
      call Message % Error(60,                                           &
                           "Missing definition of porous region "    //  &
                           trim(porous_region_rank) // ".  "         //  &
                           "Check its definition in the control file ",  &
                           file=__FILE__, line=__LINE__)
    end if

  end do  ! reg
# endif

  !---------------------------!
  !                           !
  !   Set porosity in cells   !
  !                           !
  !---------------------------!
  do reg = 1, Por % n_regions
    call Por % Set_Porosity_In_Cells(Grid, reg)
  end do

  ! Mahir's correction
  do c = 1, Grid % n_cells
    Grid % por(c) = 0
    do reg = 1, Por % n_regions
      if (Por % region(reg) % cell_porous(c)) then
        Grid % por(c) = reg
        exit
      end if
    end do
  end do

  end subroutine
