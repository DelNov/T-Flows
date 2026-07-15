!==============================================================================!
  subroutine Insert_Layers_Driver(Convert, Grid, g, n_grids)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type)     :: Convert
  type(Grid_Type), target :: Grid
  integer,     intent(in) :: g
  integer,     intent(in) :: n_grids
!-----------------------------------[Locals]-----------------------------------!
  character(SL),        save :: answer
  integer                    :: i_reg, i_lay
  integer,              save :: n_layer_regions
  integer,              save :: n_layers
  logical                    :: has_positive, has_negative
  integer, allocatable, save :: layer_regions(:)
  real,    allocatable, save :: layer_thickness(:)
!==============================================================================!

  !---------------------------------------!
  !   Print regions list, important for   !
  !    the user to see what is defined    !
  !---------------------------------------!
  call Print_Regions_List(Grid)

  !---------------------------------------------------------------!
  !   If this is the first call, ask user about boundary layers   !
  !---------------------------------------------------------------!
  if(g .eq. 1) then

    call Message % Framed(60,                                      &
      "Inserting boundary layers.                            ",    &
      "Type ordinal number(s) of boundary-condition regions. " //  &
      "\n \n                                                 " //  &
      "Example:  2  5  7                                     " //  &
      "\n \n                                                 " //  &
      "Then type boundary-layer thicknesses on the next line." //  &
      "Positive values extrude layers out of the domain. \n  " //  &
      "Negative values insert layers into the domain.        " //  &
      " \n \n                                                " //  &
      "Extruding example:  0.004  0.003  0.002  0.001 \n     " //  &
      "Intruding example: -0.004 -0.003 -0.002 -0.001        " //  &
      "\n \n                                                 " //  &
      "Type skip to skip insertion of boundary layers.")

    ! Read number of regions or skip
    call File % Read_Line(5, key_log_entry =  &
                          "# Boundary-layer regions (or skip)")
    answer = Line % tokens(1)
    call String % To_Upper_Case(answer)

    ! User wants to skip insertion of boundary layers, fine
    if(answer .eq. 'SKIP') then
      print *, "# Skipping insertion of boundary layers."
      return
    end if

    ! First line is the list of boundary regions to become boundary layers.
    ! Number of tokes corresponds to the number of regions
    n_layer_regions = Line % n_tokens
    call Enlarge % Array_Int(layer_regions, (/1, n_layer_regions/))

    ! Read all the regions
    do i_reg = 1, n_layer_regions
      read(Line % tokens(i_reg), *) layer_regions(i_reg)
      Assert(layer_regions(i_reg) .ne. 0)
    end do

    ! Read the thickness of each boundary layer
    call File % Read_Line(5, key_log_entry = "# Boundary-layer thicknesses")

    ! Number of tokes corresponds to the number of layers
    n_layers = Line % n_tokens
    call Enlarge % Array_Real(layer_thickness, (/1, n_layers/))

    ! Read all the layers
    do i_lay = 1, n_layers
      read(Line % tokens(i_lay), *) layer_thickness(i_lay)
    end do

    ! Check if thicknesses are positive or negative
    has_positive = any(layer_thickness(1:n_layers) .gt. 0.0)
    has_negative = any(layer_thickness(1:n_layers) .lt. 0.0)

    ! Make sure the sign of thicknesses is not mixed
    if(has_positive .and. has_negative) then
      call Message % Error(60,  &
        "Boundary layer thicknesses cannot have mixed signs. " //  &
        "Use positive values to extrude layers, or negative "  //  &
        "values to insert them.  This error is critical!",         &
        file=__FILE__, line=__LINE__)
    end if

  end if

  ! Set all layer thicknesses to be positive from this point on
  layer_thickness(1:n_layers) = abs(layer_thickness(1:n_layers))

  !---------------------------------------------------------------------------!
  !   If this is the first run (the first grid), and you are stil here (the   !
  !   user didn't choose to skip insertion of boundary layers, compress the   !
  !   domain for the thickness of all boundary layers.                        !
  !---------------------------------------------------------------------------!
  if(g .eq. 1 .and. has_negative) then
    call Convert % Insert_Layers(Grid,                 &
                                 n_layer_regions,      &
                                 layer_regions,        &
                                 n_layers,             &
                                 layer_thickness,      &
                                 compress_domain = .true.)
  end if

  !-----------------------------------------------------------------------!
  !   If this is the last grid, add boundary layers to selected regions   !
  !   Not that, this can be either the primal (if n_grids .eq. 1) or      !
  !   dual grid (if n_grids .eq. 2).  If dual grid is created, boundary   !
  !   layers are not added to the primal grid.                            !
  !-----------------------------------------------------------------------!
  if(g .eq. n_grids) then
    call Convert % Insert_Layers(Grid,                 &
                                 n_layer_regions,      &
                                 layer_regions,        &
                                 n_layers,             &
                                 layer_thickness,      &
                                 compress_domain = .false.)
  end if

  end subroutine
