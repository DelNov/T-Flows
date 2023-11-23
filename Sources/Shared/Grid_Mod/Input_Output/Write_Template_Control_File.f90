!==============================================================================!
  subroutine Write_Template_Control_File(Grid)
!------------------------------------------------------------------------------!
!>  The subroutine is invoked from Convert and Generate, and it creates a
!>  template control file tailored to the grid's boundary conditions for use in
!>  Process. This template file includes all mandatory components essential for
!>  Process to operate. Additionally, it contains numerous instructive comments,
!>  providing a solid foundation for users to further develop and modify.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! parent, grid for which template
                            !! control file being created
!-----------------------------------[Locals]-----------------------------------!
  integer       :: j, fu
  character(SL) :: work
  character(1)  :: CR = char(10)
!==============================================================================!

  work = Grid % name
  call String % To_Lower_Case(work)
  call File % Open_For_Writing_Ascii("control_template_for_"//trim(work), fu)

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# This is a minimum control file for T-Flows for the domain'
  write(fu,'(a)') '# "'//trim(work)//'"'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# One can use it as a staring point for defining entries'
  write(fu,'(a)') '# in T-Flows'' control file.  These files are created at'
  write(fu,'(a)') '# the end of each call to Conver or Generate sub-programs.'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# Each line starting with a # is, obviously, a comment.'
  write(fu,'(a)') '# Lines starting with ! or % are also skipped as comments.'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a1)', advance='no') CR

  work = Grid % name
  call String % To_Lower_Case(work)
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Problem name must be specified and corresponds to the'
  write(fu,'(a)') '# base name (without extensions) of the grid file used.'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a,a)') '  PROBLEM_NAME        ', trim(work)
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Boundary conditions also must be specified and are'
  write(fu,'(a)') '# listed here, for this grid, with some dummy values'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# T-Flows neglects what it doesn''t need to read, so in'
  write(fu,'(a)') '# the lines below, only "wall" will be read as the type of'
  write(fu,'(a)') '# boundary conditions, everything in () braces is skipped.'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# The same holds for the values listed.  If simulation is'
  write(fu,'(a)') '# laminar or LES, T-Flows will skip the values for "k",'
  write(fu,'(a)') '# "eps", "zeta" and "f22"'
  write(fu,'(a)') '#-----------------------------------------------------------'

  do j = Boundary_Regions()
    work = Grid % region % name(j)
    call String % To_Lower_Case(work)
    write(fu,'(a,a)') '  BOUNDARY_CONDITION ', trim(work)
    write(fu,'(a)') '    TYPE             wall  (or: '  //  &
                    'inflow / outflow / pressure / convective)'
    write(fu,'(a)') '    VARIABLES        u     v     w     t    '  //  &
                                   'kin   eps    zeta   f22'
    write(fu,'(a)') '    VALUES           0.0   0.0   0.0   10   '  //  &
                                   '0.0   1e-3   0.0    1e-3'
    write(fu,'(a1)', advance='no') CR
  end do

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# And, that''s it, what is listed above is the bare minimum.'
  write(fu,'(a)') '# If nothing else is prescribed, T-Flows will take default'
  write(fu,'(a)') '# values as documented in Documents/all_control_keywords.'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# Clearly, such a simulation wouldn''t be very useful, so'
  write(fu,'(a)') '# in the following are some of the most essential options'
  write(fu,'(a)') '# for the control file. Simply uncomment what you need.'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# All the entries (except the sub-entries such as: TYPE,'
  write(fu,'(a)') '# VARIABLES/VALUES after each of the BOUNDARY_CONDITION'
  write(fu,'(a)') '# can be inserted in any order.  T-Flows will take whatever'
  write(fu,'(a)') '# it needs, whenever it needs it.  The entire boundary'
  write(fu,'(a)') '# section may as well be at he end of the file.'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Time stepping'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '! TIME_STEP                 0.001'
  write(fu,'(a)') '! NUMBER_OF_TIME_STEPS   1200'
  write(fu,'(a)') '! RESULTS_SAVE_INTERVAL   300  (how often will it save)'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Physical properties'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '! MASS_DENSITY           1.0'
  write(fu,'(a)') '! THERMAL_CONDUCTIVITY   1.4e-4'
  write(fu,'(a)') '! DYNAMIC_VISCOSITY      1.0e-4'
  write(fu,'(a)') '! HEAT_CAPACITY          1.0'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Monitoring points'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '! NUMBER_OF_MONITORING_POINTS  4'
  write(fu,'(a)') '!   MONITORING_POINT_001       0.0  0.0  2.0'
  write(fu,'(a)') '!   MONITORING_POINT_002       0.0  0.0  2.5'
  write(fu,'(a)') '!   MONITORING_POINT_003       0.0  0.0  3.0'
  write(fu,'(a)') '!   MONITORING_POINT_004       0.0  0.0  3.5'
  write(fu,'(a)') '! POINT_FOR_MONITORING_PLANES  0.17  0.17  1.51'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# Initial condition'
  write(fu,'(a)') '#'
  write(fu,'(a)') '# These are very simular to BOUNDARY_CONDITION section above'
  write(fu,'(a)') '# and also here sub-entries VARIABLES and VALUES must be'
  write(fu,'(a)') '# specified in the given order.'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '! INITIAL_CONDITION'
  write(fu,'(a)') '!   VARIABLES        u     v     w     t    '  //  &
                                   'kin    eps    zeta     f22'
  write(fu,'(a)') '!   VALUES           0.0   0.0   0.0   10   '  //  &
                                   '1e-2   1e-3   6.6e-2   1e-3'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# For flows with inlets and outlets, initialization by a'
  write(fu,'(a)') '# pressure-like potential field could prove to be useful'
  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '!  POTENTIAL_INITIALIZATION     yes'
  write(fu,'(a1)', advance='no') CR

  write(fu,'(a)') '#-----------------------------------------------------------'
  write(fu,'(a)') '# For more customizations, check the file:'
  write(fu,'(a)') '# Documents/all_control_keywords.'
  write(fu,'(a)') '#-----------------------------------------------------------'

  close(fu)

  end subroutine
