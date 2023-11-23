!==============================================================================!
  subroutine Create_Porosity(Por, Grid)
!---------------------------------[Arguments]----------------------------------!
  class(Porosity_Type)    :: Por
  type(Grid_Type), target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: reg, n
  logical       :: found
  character(SL) :: porous_region_rank
  character(SL) :: next_strings(MAX_STRING_ITEMS)
! Muhamed doesn't need those and Bojan forgot what they were
! real          :: c1c2(2), def(2) = (/0.0, 0.0/)
!==============================================================================!

  ! Read number of Porous regions from control file
  call Control % Read_Int_Item('NUMBER_OF_POROUS_REGIONS', 0, &
                               Por % n_regions, .true.)

  if(Por % n_regions .eq. 0) return

  Por % pnt_grid => Grid

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  allocate(Por % region(Por % n_regions))

  !----------------------------------------!
  !                                        !
  !   Read porous regions characterisics   !
  !                                        !
  !----------------------------------------!
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

      ! Muhamed doesn't need those and Bojan forgot what they were
      ! call Control % Read_Real_Vector('C1X_C2X', 2, def, c1c2, .true.)
      ! Por % region(reg) % c1_x = c1c2(1)
      ! Por % region(reg) % c2_x = c1c2(2)
      ! call Control % Read_Real_Vector('C1Y_C2Y', 2, def, c1c2, .true.)
      ! Por % region(reg) % c1_y = c1c2(1)
      ! Por % region(reg) % c2_y = c1c2(2)
      ! call Control % Read_Real_Vector('C1Z_C2Z', 2, def, c1c2, .true.)
      ! Por % region(reg) % c1_z = c1c2(1)
      ! Por % region(reg) % c2_z = c1c2(2)
    else
      call Message % Error(60,                                           &
                           "Missing definition of porous region "    //  &
                           trim(porous_region_rank) // ".  "         //  &
                           "Check its definition in the control file ",  &
                           file=__FILE__, line=__LINE__)
    end if

  end do  ! reg

  !---------------------------!
  !                           !
  !   Set porosity in cells   !
  !                           !
  !---------------------------!
  do reg = 1, Por % n_regions
    call Por % Set_Porosity_In_Cells(Grid, reg)
  end do

  end subroutine
