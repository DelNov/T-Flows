!==============================================================================!
  subroutine Create_Porosity(Por, Grid)
!---------------------------------[Arguments]----------------------------------!
  class(Porosity_Type)    :: Por
  type(Grid_Type), target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: reg
  logical       :: found
  character(SL) :: porous_region_rank
  real          :: c1c2(2), def(2) = (/0.0, 0.0/)
!==============================================================================!

  ! Read number of Porous regions from control file
  call Control_Mod_Read_Int_Item('NUMBER_OF_POROUS_REGIONS', 0, &
                                  Por % n_regions, .true.)

  if(Por % n_regions .eq. 0) return

  Por % pnt_grid => Grid

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  allocate(Por % region(Por % n_regions))
  do reg = 1, Por % n_regions
    allocate(Por % region(reg) % cell_porous(Grid % n_cells))
    Por % region(reg) % cell_porous(:) = .false.
  end do

  !----------------------------------------!
  !                                        !
  !   Read porous regions characterisics   !
  !                                        !
  !----------------------------------------!
  do reg = 1, Por % n_regions

    ! Set region's name
    write(porous_region_rank, '(a,i3.3)') 'POROUS_REGION_', reg

    ! Look for it
    call Control_Mod_Position_At_One_Key(porous_region_rank, found, .true.)

    !------------------------------------------!
    !   Found the section with porous region   !
    !------------------------------------------!
    if (found) then

      ! Read "on", otherwise you will always find the first mention of STL_FILE
      call Control_Mod_Read_Char_Item_On('STL_FILE',  'default.stl',    &
                                         Por % region(reg) % stl_name,  &
                                         .true.)

      print '(a)', trim(Por % region(reg) % stl_name)

      call Control_Mod_Read_Real_Array('C1X_C2X', 2, def, c1c2, .true.)
      Por % region(reg) % c1_x = c1c2(1)
      Por % region(reg) % c2_x = c1c2(2)
      call Control_Mod_Read_Real_Array('C1Y_C2Y', 2, def, c1c2, .true.)
      Por % region(reg) % c1_y = c1c2(1)
      Por % region(reg) % c2_y = c1c2(2)
      call Control_Mod_Read_Real_Array('C1Z_C2Z', 2, def, c1c2, .true.)
      Por % region(reg) % c1_z = c1c2(1)
      Por % region(reg) % c2_z = c1c2(2)
    else
      call Message % Error(44,                                           &
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
