!==============================================================================!
  subroutine Write_Dimensions_Info(base, block)
!------------------------------------------------------------------------------!
!   Writes info on dimensional data in DB
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  integer           :: base, block ! base, block
!----------------------------------[Locals]------------------------------------!
  integer       :: base_id     ! base index number
  integer       :: block_id    ! block index number
  integer       :: coord_id    ! coordinate index number
  integer       :: solution_id ! solution index
  integer       :: n_fields    ! number of fields in solution node
  character(SL) :: fieldname   ! name of field node
  integer       :: data_type   ! RealDouble
  real          :: exponents(5)! dimensional exponents
  integer       :: c
  logical       :: dimensional_bool
  integer       :: error
!==============================================================================!

  ! Set input parameters
  base_id     = base
  block_id    = block
  solution_id = 1

  !----------------------------!
  !                            !
  !   Write DimensionalUnits   !
  !                            !
  !----------------------------!

  ! Go to CGNSBase_t DB node
  call Cg_Goto_F(file_id,  & !(in )
                 base_id,  & !(in )
                 error,    & !(out)
                 'end')      !(out)

  ! Write Units
  call Cg_Units_Write_F(Kilogram,  & !(in )
                        Meter,     & !(in )
                        Second,    & !(in )
                        Kelvin,    & !(in )
                        Degree,    & !(in )
                        error)       !(out)

  !--------------------------------!
  !                                !
  !   Write DimensionalExponents   !
  !                                !
  !--------------------------------!
  
  !-----------------!
  !   Coordinates   !
  !-----------------!

  do coord_id = 1, cgns_base(base) % cell_dim

    ! Go to GridCoordinates_t / DataArray_t DB node
    call Cg_Goto_F(file_id,              & !(in )
                   base_id,              & !(in )
                   error,                & !(out)
                   'Zone_t',             & !(in )
                   block_id,             & !(in )
                   'GridCoordinates_t',  & !(in )
                   1,                    & !(in )
                   'DataArray_t',        & !(in )
                   coord_id,             & !(in )
                   'end')                  !(in )

    if (error .ne. 0) then
      print *, '#         Failed to navigate to GridCoordinates_t subnodes'
      call Cg_Error_Exit_F()
    end if

    ! Write DataClass
    call Cg_Dataclass_Write_F(Dimensional,  & !(in )
                              error)          !(out)
    if (error .ne. 0) then
      print *, '#         Failed to write DataClass'
      call Cg_Error_Exit_F()
    end if

    ! TO DO: pull it from Work_Mod
    exponents    = 0.
    exponents(2) = 1.

    ! Write DimensionalExponents
    call Cg_Exponents_Write_F(RealDouble,  & !(in )
                              exponents,   & !(in )
                              error)         !(out)

    if (error .ne. 0) then
      print *, '#         Failed to write DimensionalExponents'
      call Cg_Error_Exit_F()
    end if

  end do

  !------------!
  !   Fields   !
  !------------!
   
  ! read n_fields
  call Cg_Nfields_F(file_id,      & !(in )
                    base_id,      & !(in )
                    block_id,     & !(in )
                    solution_id,  & !(in )
                    n_fields,     & !(out)
                    error)          !(out)

    if (error .ne. 0) then
      print *, '#         Failed to read number of fields'
      call Cg_Error_Exit_F()
    end if

  do c = 1, n_fields

    call Cg_Field_Info_F(file_id,      & !(in )
                         base_id,      & !(in )
                         block_id,     & !(in )
                         solution_id,  & !(in )
                         c,            & !(in )
                         data_type,    & !(out)
                         fieldname,    & !(out)
                         error)          !(out)

    if (error .ne. 0) then
      print *, '#         Failed to read fieldname'
      call Cg_Error_Exit_F()
    end if

    ! TO DO: pull it from Work_Mod
    exponents   =0.
    dimensional_bool = .true.

    if (fieldname .eq. 'Density') then
      exponents(1)=1.
      exponents(2)=-3.
    else if (fieldname(1:8) .eq. 'Velocity') then
      exponents(2)=1.
      exponents(3)=-1.
    else if (fieldname .eq. 'Pressure') then
      exponents(1)=1.
      exponents(2)=-1.
      exponents(3)=-2.
    else if (fieldname .eq. 'Temperature') then
      exponents(4)=1.
    else if (fieldname .eq. 'MeanTemperature') then
      exponents(4)=1.
    else if (fieldname .eq. 'MeanTemperature') then
      exponents(4)=1.
    else if (fieldname .eq. 'TemperatureFluctuations') then
      exponents(4)=2.
    else if (fieldname(1:17) .eq. 'TurbulentHeatFlux') then
      exponents(2)=1.
      exponents(3)=-1.
      exponents(4)=1.
    else if (fieldname .eq. 'TurbulentEnergyKinetic') then
      exponents(2)=2.
      exponents(3)=-2.
    else if (fieldname .eq. 'TurbulentDissipation') then
      exponents(2)=2.
      exponents(3)=-3.
    else if (fieldname .eq. 'TurbulentEnergyKineticProduction') then
      exponents(2)=2.
      exponents(3)=-3.
    else if (fieldname .eq. 'TurbulentQuantityV2') then
      exponents(2)=2.
      exponents(3)=-2.
    else if (fieldname .eq. 'TurbulentViscosity') then
      exponents(2)=2.
      exponents(3)=-1.
    else if (fieldname .eq. 'VorticityMagnitude') then
      exponents(3)=-1.
    else if (fieldname(1:14) .eq. 'ReynoldsStress') then
      exponents(2)=2.
      exponents(3)=-2.
    else if (fieldname .eq. 'WallDistance') then
      exponents(2)=1.
    else ! all others are NonDimensional
      dimensional_bool = .false.
    end if

    ! Go to FlowSolution_t / DataArray_t DB node
    call Cg_Goto_F(file_id,           & !(in )
                   base_id,           & !(in )
                   error,             & !(out)
                   'Zone_t',          & !(in )
                   block_id,          & !(in )
                   'FlowSolution_t',  & !(in )
                   1,                 & !(in )
                   'DataArray_t',     & !(in )
                   c,                 & !(in )
                   'end')

    if (error .ne. 0) then
      print *, '#         Failed to navigate to FlowSolution_t subnodes'
      call Cg_Error_Exit_F()
    end if

    if (dimensional_bool) then

      ! Write DataClass
      call Cg_Dataclass_Write_F(Dimensional,  & !(in )
                                error)          !(out)
      if (error .ne. 0) then
        print *, '#         Failed to write DataClass'
        call Cg_Error_Exit_F()
      end if

      ! Write DimensionalExponents
      call Cg_Exponents_Write_F(RealDouble,  & !(in )
                                exponents,   & !(in )
                                error)         !(out)
      if (error .ne. 0) then
        print *, '#         Failed to write DimensionalExponents'
        call Cg_Error_Exit_F()
      end if

    else

      ! Write DataClass
      call Cg_Dataclass_Write_F(NondimensionalParameter,  & !(in )
                                error)                      !(out)
      if (error .ne. 0) then
        print *, '#         Failed to write DataClass'
        call Cg_Error_Exit_F()
      end if

    end if

  end do

  end subroutine
