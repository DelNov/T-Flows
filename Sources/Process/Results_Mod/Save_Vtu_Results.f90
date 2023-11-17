!==============================================================================!
  subroutine Save_Vtu_Results(Results, Flow, Turb, Vof, Swarm,  &
                              plot_inside, domain)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)         :: Results
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  logical                     :: plot_inside  ! plot results inside?
  integer,           optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Var_Type),  pointer     :: phi
  integer(SP)                  :: data_size
  integer                      :: data_offset, cell_offset, i_fac, reg
  integer                      :: s, n, n_conns, n_polyg, sc, f8, f9, run
  integer                      :: s1, s2, c1, c2, c_f, c_l
  real                         :: dist1, dist2
  character(SL)                :: name_out_8, name_out_9, name_mean
  character(SL)                :: str1, str2, str_time, str_var
  integer, pointer, contiguous :: int_save(:), type_save(:), offs_save(:)
  real,    pointer, contiguous :: save_01(:), save_02(:), save_03(:)
  real,    pointer, contiguous :: save_04(:), save_05(:), save_06(:)
  real,    pointer, contiguous :: var_ins(:)
  real,    pointer, contiguous :: v2_calc(:), kin_vis_t(:), phi_save(:)
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

  call Profiler % Start('Save_Vtu_Results')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Take aliases
  Grid => Flow % pnt_grid

  if(.not. plot_inside .and. .not. Results % boundary) return

  call Work % Connect_Int_Cell(int_save, type_save, offs_save)
  call Work % Connect_Real_Cell(save_01, save_02, save_03,  &
                                save_04, save_05, save_06)
  call Work % Connect_Real_Cell(var_ins, v2_calc, kin_vis_t, phi_save)

  !------------------------------------------------!
  !   Mark the beginnings and end of cell ranges   !
  !------------------------------------------------!

  ! For cells inside
  if(plot_inside) then
    c_f = 1
    c_l = Grid % n_cells
    if(.not. PLOT_BUFFERS) c_l = Grid % n_cells - Grid % Comm % n_buff_cells

  ! For boundary cells
  else
    c_f = -Grid % n_bnd_cells
    c_l = -1
    if(.not. PLOT_BUFFERS) then
      do c_f = -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_f) .eq. This_Proc()) exit
      end do
      do c_l = -1, -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_l) .eq. This_Proc()) exit
      end do
    end if
  end if

  !-------------------------------------------------------------------------!
  !   Count connections and polygons in this Grid, you will need it later   !
  !-------------------------------------------------------------------------!
  n_conns = 0
  n_polyg = 0
  if(plot_inside) then
    ! Connections
    do c1 = c_f, c_l
      n_conns = n_conns + abs(Grid % cells_n_nodes(c1))
    end do
    ! Polygons
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                   ! add one for number of polyfs
        do i_fac = 1, Grid % cells_n_faces(c1)  ! add all faces and their nodes
          s = Grid % cells_f(i_fac, c1)
          n = Grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end do
  else
    do c2 = c_f, c_l
      n_conns = n_conns + Grid % cells_n_nodes(c2)
    end do
  end if

  call Global % Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File % Set_Name(name_out_8,                    &
                         time_step = Time % Curr_Dt(),  &
                         extension = '.pvtu',           &
                         domain    = domain)
    call File % Set_Name(name_out_9,                              &
                         processor = (/This_Proc(), N_Procs()/),  &
                         time_step = Time % Curr_Dt(),            &
                         extension = '.vtu',                      &
                         domain    = domain)
  else
    call File % Set_Name(name_out_8,                    &
                         time_step = Time % Curr_Dt(),  &
                         appendix  = '-bnd',            &
                         extension = '.pvtu',           &
                         domain=domain)
    call File % Set_Name(name_out_9,                              &
                         processor = (/This_Proc(), N_Procs()/),  &
                         time_step = Time % Curr_Dt(),            &
                         appendix  = '-bnd',                      &
                         extension = '.vtu',                      &
                         domain    = domain)
  end if

  if(Parallel_Run() .and. First_Proc()) then
    call File % Open_For_Writing_Binary(name_out_8, f8)
  end if
  call File % Open_For_Writing_Binary(name_out_9, f9)

  write(str_time,'(E16.9)')Time % Get_Time()
  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(Parallel_Run() .and. First_Proc())  then
    write(f8) IN_0 // '<?xml version="1.0"?>'              // LF
    write(f8) IN_0 // '<VTKFile type="PUnstructuredGrid">' // LF
    write(f8) IN_1 // '<PUnstructuredGrid GhostLevel="1">' // LF
    write(f8) IN_2 // '<PFieldData>'                       // LF
    ! TIME must be capitalized for visit
    write(f8) IN_3 // '<PDataArray type="Float64" Name="TIME" ' // &
                  'NumberOfTuples="1" format="ascii"> ' // trim(str_time) // LF
    write(f8) IN_3 // '</PDataArray>'                      // LF
    write(f8) IN_2 // '</PFieldData>'                      // LF
  end if

  write(f9) IN_0 // '<?xml version="1.0"?>'                           // LF
  write(f9) IN_0 // '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                    'byte_order="LittleEndian">'                      // LF
  write(f9) IN_1 // '<UnstructuredGrid>'                              // LF
  write(f9) IN_2 // '<FieldData>'                                     // LF
  ! TIME must be capitalized for visit
  write(f9) IN_3 // '<DataArray type="Float64" Name="TIME" ' // & 
                    'NumberOfTuples="1" format="ascii">' // trim(str_time) // LF
  write(f9) IN_3 // '</DataArray>'                                    // LF
  write(f9) IN_2 // '</FieldData>'                                    // LF

  write(str1,'(i0.0)') Grid % n_nodes
  if(plot_inside) then
    write(str2,'(i0.0)') (c_l-c_f+1)
  else
    if((c_l-c_f+1) .eq. 0) then
      write(str2,'(i1)')   (c_l-c_f+1)  ! 0.0 doesn't work for zero :-/
    else
      write(str2,'(i0.0)') (c_l-c_f+1)
    end if
  end if
  write(f9) IN_2 // '<Piece NumberOfPoints="' // trim(str1)    //  &
                    '" NumberOfCells ="' // trim(str2) // '">' // LF

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  if(Parallel_Run() .and. First_Proc())  then
    write(f8) IN_3 // '<PPoints>' // LF
    write(f8) IN_4 // '<PDataArray type='//floatp  //  &
                      ' NumberOfComponents="3"/>'  // LF
    write(f8) IN_3 // '</PPoints>' // LF
  end if

  write(str1, '(i1)') data_offset
  write(f9) IN_3 // '<Points>'                        // LF
  write(f9) IN_4 // '<DataArray type='//floatp        //  &
                    ' NumberOfComponents="3"'         //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
  write(f9) IN_4 // '</DataArray>' // LF
  write(f9) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3

  !-----------!
  !   Cells   !
  !-----------!
  write(f9) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(f9) IN_4 // '<DataArray type='//intp          //  &
                    ' Name="connectivity"'            //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
  write(f9) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Fill up an array with cell offsets and save the header only
  cell_offset = 0
  if(plot_inside) then
    do c1 = c_f, c_l
      cell_offset   = cell_offset + abs(Grid % cells_n_nodes(c1))
      offs_save(c1) = cell_offset
    end do
  else
    do c2 = c_f, c_l
      cell_offset   = cell_offset + Grid % cells_n_nodes(c2)
      offs_save(c2) = cell_offset
    end do
  end if
  call Results % Save_Vtu_Scalar_Int("offsets", plot_inside,  &
                                     offs_save(c_f:c_l),      &
                                     f8, f9, data_offset, 1)  ! 1 => header only

  ! Fill up an array with cell types and save the header only
  if(plot_inside) then
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .eq. 8) type_save(c1) = VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c1) .eq. 6) type_save(c1) = VTK_WEDGE
      if(Grid % cells_n_nodes(c1) .eq. 4) type_save(c1) = VTK_TETRA
      if(Grid % cells_n_nodes(c1) .eq. 5) type_save(c1) = VTK_PYRAMID
      if(Grid % cells_n_nodes(c1) .lt. 0) type_save(c1) = VTK_POLYHEDRON
    end do
  else
    do c2 = c_f, c_l
      if(Grid % cells_n_nodes(c2) .eq. 4) then
        type_save(c2) = VTK_QUAD
      else if(Grid % cells_n_nodes(c2) .eq. 3) then
        type_save(c2) = VTK_TRIANGLE
      else
        type_save(c2) = VTK_POLYGON
      end if
    end do
  end if
  call Results % Save_Vtu_Scalar_Int("types", plot_inside,  &
                                     type_save(c_f:c_l),    &
                                     f8, f9, data_offset, 1)  ! 1 => header only

  ! Write parts of header for polyhedral cells
  if(n_polyg > 0) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next


    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faceoffsets"'            //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + (c_l-c_f+1) * IP  ! prepare for next

  end if

  write(f9) IN_3 // '</Cells>'     // LF

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  if(Parallel_Run() .and. First_Proc())  then
    write(f8) IN_3 // '<PCellData>' // LF
  end if
  write(f9) IN_3 // '<CellData>' // LF

  !----------------!
  !                !
  !   Two sweeps   !
  !                !
  !----------------!
  do run = 1, 2

    !------------------------------------------!
    !   Save remnants of the Grid definition   !
    !------------------------------------------!
    if(run .eq. 2) then

      ! Save the nodes' coordinates
      data_size = int(Grid % n_nodes * RP * 3, SP)
      write(f9) data_size
      do n = 1, Grid % n_nodes
        write(f9) Grid % xn(n), Grid % yn(n), Grid % zn(n)
      end do

      ! Save connections
      data_size = int(n_conns * IP, SP)
      write(f9) data_size
      if(plot_inside) then
        do c1 = c_f, c_l

          ! Tetrahedral, pyramid, wedge and hexahedral cells
          if( any( Grid % cells_n_nodes(c1) .eq. (/4,5,6,8/)  ) ) then
            write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c1), c1))-1

          ! Polyhedral cells
          else if(Grid % cells_n_nodes(c1) .lt. 0) then
            write(f9) (Grid % cells_n(1:-Grid % cells_n_nodes(c1), c1))-1

          end if
        end do
      else  ! plot only boundary
        do c2 = c_f, c_l

          ! All cell types
          write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c2), c2))-1

        end do
      end if

      ! Save cell offsets
      call Results % Save_Vtu_Scalar_Int("offsets", plot_inside,    &
                                         offs_save(c_f:c_l),        &
                                         f8, f9, data_offset, run)
      ! Save cell types
      call Results % Save_Vtu_Scalar_Int("types", plot_inside,      &
                                         type_save(c_f:c_l),        &
                                         f8, f9, data_offset, run)

      if(n_polyg > 0) then

        ! Write polyhedral cells' faces
        data_size = int(n_polyg * IP, SP)
        write(f9) data_size

        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            write(f9) Grid % cells_n_faces(c1)      ! write number of polygons
            do i_fac = 1, Grid % cells_n_faces(c1)  ! write nodes of each polygn
              s = Grid % cells_f(i_fac, c1)
              if(Grid % faces_s(s) .ne. 0) then  ! face has a shadow, if it ...
                s1 = s                           ! ... is closer, plot that!
                s2 = Grid % faces_s(s)
                dist1 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s1), Grid % yf(s1), Grid % zf(s1))
                dist2 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
                if(dist1 < dist2) s = s1
                if(dist2 < dist1) s = s2
              end if
              n = Grid % faces_n_nodes(s)
              write(f9) n, (Grid % faces_n(1:n, s))-1
            end do
          end if
        end do

        ! Write polyhedral cells' faces offsets
        data_size = int((c_l-c_f+1) * IP, SP)
        write(f9) data_size

        cell_offset = 0
        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            cell_offset = cell_offset + 1           ! to store number of polygs
            do i_fac = 1, Grid % cells_n_faces(c1)  ! to store all the nodes
              s = Grid % cells_f(i_fac, c1)         ! of each polygon
              n = Grid % faces_n_nodes(s)
              cell_offset = cell_offset + 1 + n
            end do
            write(f9) cell_offset

          ! Not a polyhedron, offsets are not needed and set to -1
          else
            write(f9) -1
          end if

        end do

      end if  ! n_polyg > 0

    end if

    !--------------------!
    !   Cell processor   !
    !--------------------!
    do c1 = c_f, c_l
      int_save(c1) = Grid % Comm % cell_proc(c1)
    end do
    do c2 = c_f, c_l
      int_save(c2) = Grid % Comm % cell_proc(c2)
    end do

    str_var =  Var_Name("Grid Processor","[1]",Results % units)
    call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                    int_save(c_f:c_l),               &
                                    f8, f9, data_offset, run)

    !-----------------!
    !   Cell thread   !
    !-----------------!
    int_save(c_f:c_l) = Grid % Vect % cell_thread(c_f:c_l)
    str_var =  Var_Name("Grid Thread","[1]",Results % units)
    call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                       int_save(c_f:c_l),            &
                                       f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(c_f:c_l) = domain
      str_var =  Var_Name("Grid Domain","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                         int_save(c_f:c_l),            &
                                         f8, f9, data_offset, run)
    end if

    !----------------------!
    !   Boundary regions   !
    !----------------------!
    if(.not. plot_inside) then
      int_save(c_f:c_l) = Grid % region % at_cell(c_f:c_l)
      str_var =  Var_Name("Boundary Condition","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Int(trim(str_var),plot_inside,    &
                                         int_save(c_f:c_l),            &
                                         f8, f9, data_offset, run)
    end if

    !--------------!
    !   Velocity   !
    !--------------!
    str_var =  Var_Name("Velocity","[m/s]",Results % units)
    call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                        Flow % u % n(c_f:c_l),        &
                                        Flow % v % n(c_f:c_l),        &
                                        Flow % w % n(c_f:c_l),        &
                                        f8, f9, data_offset, run)
    !--------------------!
    !   Courant number   !
    !--------------------!
    if(plot_inside) then
      call Flow % Calculate_Courant_In_Cells(save_01)
      str_var =  Var_Name("Courant Number","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                          save_01(c_f:c_l),            &
                                          f8, f9, data_offset, run)
    end if

    !---------------!
    !   Potential   !
    !---------------!
    str_var =  Var_Name("Potential","[m^2/s]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Flow % potential(c_f:c_l),    &
                                        f8, f9, data_offset, run)

    !--------------------------------------!
    !   Pressure correction and pressure   !
    !--------------------------------------!
    str_var =  Var_Name("Pressure Correction","[Pa]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),            &
                                        plot_inside,              &
                                        Flow % pp % n(c_f:c_l),   &
                                        f8, f9, data_offset, run)
    save_01(:) = 0.0
    save_02(:) = 0.0
    save_03(:) = 0.0
    do c1 = c_f, c_l
      save_01(c1) = Flow % pp % x(c1) * Grid % vol(c1)
      save_02(c1) = Flow % pp % y(c1) * Grid % vol(c1)
      save_03(c1) = Flow % pp % z(c1) * Grid % vol(c1)
    end do
    str_var =  Var_Name("Pressure Correction Force","[N]",Results % units)
    call Results % Save_Vtu_Vector_Real(trim(str_var),      &
                                        plot_inside,        &
                                        save_01(c_f:c_l),   &
                                        save_02(c_f:c_l),   &
                                        save_03(c_f:c_l),   &
                                        f8, f9, data_offset, run)

    str_var =  Var_Name("Pressure","[Pa]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                        Flow % p % n(c_f:c_l),         &
                                        f8, f9, data_offset, run)
    save_01(:) = 0.0
    save_02(:) = 0.0
    save_03(:) = 0.0
    do c1 = c_f, c_l
      save_01(c1) = Flow % p % x(c1) * Grid % vol(c1)
      save_02(c1) = Flow % p % y(c1) * Grid % vol(c1)
      save_03(c1) = Flow % p % z(c1) * Grid % vol(c1)
    end do

    str_var =  Var_Name("PressureForce","[N]",Results % units)
    call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                        save_01(c_f:c_l),             &
                                        save_02(c_f:c_l),             &
                                        save_03(c_f:c_l),             &
                                        f8, f9, data_offset, run)

    !-----------------!
    !   Temperature   !
    !-----------------!
    if(Flow % heat_transfer) then
      str_var =  Var_Name("Temperature","[K]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Flow % t % n(c_f:c_l),        &
                                          f8, f9, data_offset, run)
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0
      do c1 = c_f, c_l
        save_01(c1) = Flow % t % x(c1)
        save_02(c1) = Flow % t % y(c1)
        save_03(c1) = Flow % t % z(c1)
      end do

      if(.not. Flow % mass_transfer) then
        call Flow % Grad_Variable(Flow % t)
      else
        call Vof % Calculate_Grad_Matrix_With_Front()
        call Vof % Grad_Variable_With_Front(Flow % t, Vof % t_sat)
        call Flow % Calculate_Grad_Matrix()
      end if

      str_var =  Var_Name("Temperature Gradients","[K/m]",Results % units)
      call Results % Save_Vtu_Vector_Real(trim(str_var),      &
                                          plot_inside,        &
                                          save_01(c_f:c_l),   &
                                          save_02(c_f:c_l),   &
                                          save_03(c_f:c_l),   &
                                          f8, f9, data_offset, run)

    end if

    !-------------------------!
    !   Physical properties   !
    !-------------------------!
    str_var =  Var_Name("Physical Density","[kg/m^3]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % density(c_f:c_l),        &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Physical Viscosity","[Pa s]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % viscosity(c_f:c_l),      &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Physical Conductivity","[W/m/K]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % conductivity(c_f:c_l),   &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Physical Capacity","[J/kg/K]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % capacity(c_f:c_l),       &
                                        f8, f9, data_offset, run)

    if(Turb % rough_walls) then
      str_var =  Var_Name("Roughness Coefficient z_o","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Turb % z_o(c_f:c_l),           &
                                          f8, f9, data_offset, run)

    end if

    !---------------------!
    !   Volume fraction   !
    !---------------------!
    if(Flow % with_interface) then
      str_var =  Var_Name("Vof Sharp","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Vof % fun % n(c_f:c_l),        &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Vof Smooth","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Vof % smooth % n(c_f:c_l),     &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Vof Curvature","[1/m]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Vof % curv(c_f:c_l),           &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Vof SurfaceNormals","[1]",Results % units)
      call Results % Save_Vtu_Vector_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Vof % nx(c_f:c_l),             &
                                          Vof % ny(c_f:c_l),             &
                                          Vof % nz(c_f:c_l),             &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Vof SurfaceTensionForce","[N]",Results % units)
      call Results % Save_Vtu_Vector_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Vof % surf_fx(c_f:c_l),        &
                                          Vof % surf_fy(c_f:c_l),        &
                                          Vof % surf_fz(c_f:c_l),        &
                                          f8, f9, data_offset, run)
      if (allocated(Vof % m_dot)) then
        str_var =  Var_Name("Vof MassTransfer","[kg/m^3/s]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var),               &
                                            plot_inside,                 &
                                            Vof % m_dot(c_f:c_l),        &
                                            f8, f9, data_offset, run)
      end if
    end if

    !---------------------------------------!
    !   Number of impacts and reflections   !
    !---------------------------------------!
    if(Flow % with_particles .and. .not. plot_inside) then
      str_var =  Var_Name("Particles Reflected","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Swarm % n_reflected(c_f:c_l),  &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Particles Deposited","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                 &
                                          plot_inside,                   &
                                          Swarm % n_deposited(c_f:c_l),  &
                                          f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      call Results % Save_Vtu_Scalar_Real(phi % name, plot_inside,   &
                                          phi % n(c_f:c_l),          &
                                          f8, f9, data_offset, run)
    end do

    !-----------------!
    !   Q-criterion   !
    !-----------------!
    phi_save(:) = 0.0
    do c1 = c_f, c_l
      phi_save(c1) = (Flow % vort(c1)**2 - Flow % shear(c1)**2)/4.
    end do
    str_var =  Var_Name("Q Criterion","[1/s^2]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),plot_inside,   &
                                        phi_save(c_f:c_l),           &
                                        f8, f9, data_offset, run)

    !--------------------------!
    !   Turbulent quantities   !
    !--------------------------!

    ! Save kin and eps
    if(Turb % model .eq. K_EPS                 .or.  &
       Turb % model .eq. K_EPS_ZETA_F          .or.  &
       Turb % model .eq. HYBRID_LES_RANS       .or.  &
       Turb % model .eq. RSM_MANCEAU_HANJALIC  .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC  ) then
      str_var =  Var_Name("Turbulent Kinetic Energy","[m^2/s^2]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                            Turb % kin % n(c_f:c_l),                    &
                            f8, f9, data_offset, run)
      str_var =  Var_Name("Turbulent Dissipation","[m^2/s^3]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                            Turb % eps % n(c_f:c_l),                    &
                            f8, f9, data_offset, run)
      str_var =  Var_Name("Turbulent Kinetic Energy Production","[m^2/s^3]",  &
                          Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                            Turb % p_kin(c_f:c_l),                      &
                            f8, f9, data_offset, run)
    end if

    ! Save zeta and f22
    if(Turb % model .eq. K_EPS_ZETA_F .or.  &
       Turb % model .eq. HYBRID_LES_RANS) then
      v2_calc(:) = 0.0
      do c1 = c_f, c_l
        v2_calc(c1) = Turb % kin % n(c1) * Turb % zeta % n(c1)
      end do
!     call Results % Save_Vtu_Scalar_Real("Turbulent Quantity V2 [m^2/s^2]",  &
!                                         plot_inside,                        &
!                                         v2_calc (c_f:c_l),                  &
!                                         f8, f9, data_offset, run)
      str_var =  Var_Name("Turbulent Quantity Zeta","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                          Turb % zeta % n(c_f:c_l),    &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Turbulent Quantity F22","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                          Turb % f22  % n(c_f:c_l),    &
                                          f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        str_var =  Var_Name("Turbulent Quantity T2","[K^2]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                            Turb % t2 % n(c_f:c_l),      &
                                            f8, f9, data_offset, run)
!       call Results % Save_Vtu_Scalar_Real("Turbulent T2 Production [K^2/s]",&
!                                           plot_inside,                      &
!                                           Turb % p_t2(c_f:c_l),             &
!                                           f8, f9, data_offset, run)
        ! unit wrong ???
        str_var =  Var_Name("Turbulent Heat Flux","[K m/s]",Results % units)
        call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,  &
                                            Turb % ut % n(c_f:c_l),      &
                                            Turb % vt % n(c_f:c_l),      &
                                            Turb % wt % n(c_f:c_l),      &
                                            f8, f9, data_offset, run)
!        call Results % Save_Vtu_Scalar_Real("Turbulent Quantity Alpha L",    &
!                                            plot_inside,                     &
!                                            Turb % alpha_l(c_f:c_l),         &
!                                            f8, f9, data_offset, run)
!        call Results % Save_Vtu_Scalar_Real("Turbulent Quantity Alpha U",    &
!                                            plot_inside,                     &
!                                            Turb % alpha_u(c_f:c_l),         &
!                                            f8, f9, data_offset, run)
      end if
    end if

    if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
      str_var =  Var_Name("Turbulent Quantity F22","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                          Turb % f22 % n(c_f:c_l),     &
                                          f8, f9, data_offset, run)
    end if

    ! Save vis and vis_t
    if(Turb % model .eq. DES_SPALART .or.  &
       Turb % model .eq. SPALART_ALLMARAS) then
      str_var =  Var_Name("Turbulent Viscosity","[Pa s]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Turb % vis % n(c_f:c_l),      &
                                          f8, f9, data_offset, run)
      str_var =  Var_Name("Vorticity Magnitude","[1/s]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Flow % vort(c_f:c_l),         &
                                          f8, f9, data_offset, run)
    end if

    kin_vis_t(:) = 0.0
    if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       Turb % model .ne. HYBRID_LES_RANS     .and.  &
       Turb % model .ne. DNS) then
      kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      str_var =  Var_Name("Eddy Over Molecular Viscosity","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                  kin_vis_t(c_f:c_l),                   &
                                  f8, f9, data_offset, run)
    end if

    if(Turb % model .eq. HYBRID_LES_RANS) then
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      str_var =  Var_Name("Rans Eddy Over Molecular Viscosity","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                  kin_vis_t(c_f:c_l),                   &
                                  f8, f9, data_offset, run)
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = Turb % vis_t_sgs(c_f:c_l) / Flow % viscosity(c_f:c_l)
      str_var =  Var_Name("Sgs Eddy Over Molecular Viscosity","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                  kin_vis_t(c_f:c_l),                   &
                                  f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

      ! Note: follows the order in which Paraview stores tensors
      str_var =  Var_Name("Reynolds Stress","[m^2/s^2]",Results % units)
      call Results % Save_Vtu_Tensor_6_Real(trim(str_var), plot_inside,   &
                                            Turb % uu % n(c_f:c_l),       &
                                            Turb % vv % n(c_f:c_l),       &
                                            Turb % ww % n(c_f:c_l),       &
                                            Turb % uv % n(c_f:c_l),       &
                                            Turb % vw % n(c_f:c_l),       &
                                            Turb % uw % n(c_f:c_l),       &
                                            f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        str_var =  Var_Name("Turbulent Heat Flux","[K m/s]",Results % units)
        call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                            Turb % ut % n(c_f:c_l),       &
                                            Turb % vt % n(c_f:c_l),       &
                                            Turb % wt % n(c_f:c_l),       &
                                            f8, f9, data_offset, run)
      end if
    end if

    ! Statistics for large-scale simulations of turbulence
    if(Turb % statistics) then
      str_var =  Var_Name("Mean Velocity","[m/s]",Results % units)
      call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                          Turb % u_mean(c_f:c_l),       &
                                          Turb % v_mean(c_f:c_l),       &
                                          Turb % w_mean(c_f:c_l),       &
                                          f8, f9, data_offset, run)
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0
      save_04(:) = 0.0
      save_05(:) = 0.0
      save_06(:) = 0.0

      ! Note: follows the order in which Paraview stores tensors
      do c1 = c_f, c_l
        save_01(c1) = Turb % uu_res(c1) - Turb % u_mean(c1) * Turb % u_mean(c1)
        save_02(c1) = Turb % vv_res(c1) - Turb % v_mean(c1) * Turb % v_mean(c1)
        save_03(c1) = Turb % ww_res(c1) - Turb % w_mean(c1) * Turb % w_mean(c1)
        save_04(c1) = Turb % uv_res(c1) - Turb % u_mean(c1) * Turb % v_mean(c1)
        save_05(c1) = Turb % vw_res(c1) - Turb % v_mean(c1) * Turb % w_mean(c1)
        save_06(c1) = Turb % uw_res(c1) - Turb % u_mean(c1) * Turb % w_mean(c1)
      end do
      str_var =  Var_Name("Mean Reynolds Stress","[m^s/s^2]",Results % units)
      call Results % Save_Vtu_Tensor_6_Real(trim(str_var), plot_inside,   &
                                        save_01(c_f:c_l),                 &
                                        save_02(c_f:c_l),                 &
                                        save_03(c_f:c_l),                 &
                                        save_04(c_f:c_l),                 &
                                        save_05(c_f:c_l),                 &
                                        save_06(c_f:c_l),                 &
                                        f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        str_var =  Var_Name("Mean Temperature","[K]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Turb % t_mean(c_f:c_l),            &
                                        f8, f9, data_offset, run)
        phi_save(:) = 0.0
        save_01(:) = 0.0
        save_02(:) = 0.0
        save_03(:) = 0.0
        do c1 = c_f, c_l
          phi_save(c1) = Turb % t2_res(c1) - Turb % t_mean(c1)*Turb % t_mean(c1)
          save_01(c1)  = Turb % ut_res(c1) - Turb % u_mean(c1)*Turb % t_mean(c1)
          save_02(c1)  = Turb % vt_res(c1) - Turb % v_mean(c1)*Turb % t_mean(c1)
          save_03(c1)  = Turb % wt_res(c1) - Turb % w_mean(c1)*Turb % t_mean(c1)
        end do
        str_var =  Var_Name("Mean Turbulent Quantity T2","[K^2]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        phi_save(c_f:c_l),                &
                                        f8, f9, data_offset, run)
        str_var =  Var_Name("Mean Turbulent Heat Flux","[K m/s]",Results % units)
        call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                        save_01(c_f:c_l),                 &
                                        save_02(c_f:c_l),                 &
                                        save_03(c_f:c_l),                 &
                                        f8, f9, data_offset, run)
      end if

      ! Scalars
      do sc = 1, Flow % n_scalars
        phi => Flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        do c1 = c_f, c_l
          phi_save(c1) = Turb % scalar_mean(sc, c1)
        end do
        call Results % Save_Vtu_Scalar_Real(name_mean, plot_inside,    &
                                            phi_save(c_f:c_l),         &
                                            f8, f9, data_offset, run)
      end do
    end if

    ! Save y+ for all turbulence models
    if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       Turb % model .ne. DNS) then
      str_var =  Var_Name("Turbulent Quantity Y Plus","[1]",Results % units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Turb % y_plus(c_f:c_l),       &
                                          f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    str_var =  Var_Name("Grid Cell Volume","[m^3]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Grid % vol(c_f:c_l),          &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Grid Cell Inertia","[m^2]",Results % units)
    call Results % Save_Vtu_Tensor_6_Real(trim(str_var), plot_inside,   &
                                          Grid % ixx(c_f:c_l),          &
                                          Grid % iyy(c_f:c_l),          &
                                          Grid % izz(c_f:c_l),          &
                                          Grid % ixy(c_f:c_l),          &
                                          Grid % iyz(c_f:c_l),          &
                                          Grid % ixz(c_f:c_l),          &
                                          f8, f9, data_offset, run)
    str_var =  Var_Name("Grid Wall Distance","[m]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Grid % wall_dist(c_f:c_l),    &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Grid Cell Delta Max","[m]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Turb % h_max(c_f:c_l),        &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Grid Cell Delta Min","[m]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Turb % h_min(c_f:c_l),        &
                                        f8, f9, data_offset, run)
    str_var =  Var_Name("Grid Cell Delta Wall","[m]",Results % units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                        Turb % h_w  (c_f:c_l),        &
                                        f8, f9, data_offset, run)

    !---------------------------------------------------------------------!
    !   Variables in the first computational point, plotted at boundary   !
    !---------------------------------------------------------------------!

    ! Engage only for boundary plots (not inside means on the boundary)
    if( .not. plot_inside ) then 

      ! Initialize working variables to zero
      save_01(:) = 0.0
      save_02(:) = 0.0
      save_03(:) = 0.0

      ! Copy internal values to boundary
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then 
          save_01(c2) = Flow % u % n(c1)
          save_02(c2) = Flow % v % n(c1)
          save_03(c2) = Flow % w % n(c1)
        end if
      end do
      call Results % Save_Vtu_Vector_Real("Velocity Near Wall [m/s]",  &
                                      plot_inside,                     &
                                      save_01(c_f:c_l),                &
                                      save_02(c_f:c_l),                &
                                      save_03(c_f:c_l),                &
                                      f8, f9, data_offset, run)

      if(Turb % model .eq. K_EPS                 .or.  &
         Turb % model .eq. K_EPS_ZETA_F          .or.  &
         Turb % model .eq. HYBRID_LES_RANS) then

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do reg = Boundary_Regions()
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            var_ins(c2) = Turb % kin % n(c1)
          end do
        end do

        str_var =  Var_Name("T.K.E. Near Wall","[m^2/s^2]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                            var_ins(c_f:c_l),             &
                                            f8, f9, data_offset, run)

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do reg = Boundary_Regions()
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            var_ins(c2) = Turb % y_plus(c1)
          end do
        end do

        str_var =  Var_Name("y+ Near Wall","[1]",Results % units)
        call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                            var_ins(c_f:c_l),           &
                                            f8, f9, data_offset, run)

        do sc = 1, Flow % n_scalars
          phi => Flow % scalar(sc)
          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do reg = Boundary_Regions()
            do s = Faces_In_Region(reg)
              c1 = Grid % faces_c(1,s)
              c2 = Grid % faces_c(2,s)

              var_ins(c2) = phi % n(c1)
            end do
          end do

          call Results % Save_Vtu_Scalar_Real("Scalar Near Wall",        &
                                              plot_inside,               &
                                              var_ins(c_f:c_l),          &
                                              f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do reg = Boundary_Regions()
            do s = Faces_In_Region(reg)
              c1 = Grid % faces_c(1,s)
              c2 = Grid % faces_c(2,s)

              var_ins(c2) = phi_save(c1)  ! Turb % scalar_mean(sc, c1)
            end do
          end do

          call Results % Save_Vtu_Scalar_Real("Mean Scalar Near Wall",  &
                                              plot_inside,              &
                                              var_ins(c_f:c_l),         &
                                              f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do reg = Boundary_Regions()
            do s = Faces_In_Region(reg)
              c1 = Grid % faces_c(1,s)
              c2 = Grid % faces_c(2,s)

              var_ins(c2) = phi % q(c2)  ! Turb % scalar_mean(sc, c1)
            end do
          end do

          call Results % Save_Vtu_Scalar_Real("Wall Scalar Flux",        &
                                              plot_inside,               &
                                              var_ins(c_f:c_l),          &
                                              f8, f9, data_offset, run)

        end do
      end if
    end if

    !----------------------!
    !                      !
    !   End of cell data   !
    !                      !
    !----------------------!
    if(run .eq. 1) then
      if(Parallel_Run() .and. First_Proc()) then
        write(f8) IN_3 // '</PCellData>' // LF
      end if
      write(f9) IN_3 // '</CellData>' // LF

      write(f9) IN_2 // '</Piece>'            // LF
      write(f9) IN_1 // '</UnstructuredGrid>' // LF

      !-------------------!
      !                   !
      !   Appended data   !
      !                   !
      !-------------------!
      write(f9) IN_0 // '<AppendedData encoding="raw">' // LF
      write(f9) '_'
    else
      write(f9) LF // IN_0 // '</AppendedData>' // LF
    end if

  end do  ! run

  write(f9) IN_0 // '</VTKFile>'          // LF
  close(f9)

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  if(Parallel_Run() .and. First_Proc()) then
    do n = 1, N_Procs()
      if(plot_inside) then
        call File % Set_Name(name_out_9,                    &
                             !processor = n,                 &
                             processor = (/n, N_Procs()/),  &
                             time_step = Time % Curr_Dt(),  &
                             extension = '.vtu',            &
                             domain    = domain)
      else
        call File % Set_Name(name_out_9,                    &
                             !/processor = n,                 &
                             processor = (/n, N_Procs()/),  &
                             time_step = Time % Curr_Dt(),  &
                             appendix  = '-bnd',            &
                             extension = '.vtu',            &
                             domain    =  domain)
      end if
      write(f8) IN_2 // '<Piece Source="', trim(name_out_9), '"/>' // LF
    end do
    write(f8) IN_1 // '</PUnstructuredGrid>' // LF
    write(f8) IN_0 // '</VTKFile>'           // LF
    close(f8)
  end if

  call Work % Disconnect_Int_Cell(int_save, type_save, offs_save)
  call Work % Disconnect_Real_Cell(save_01, save_02, save_03,  &
                                   save_04, save_05, save_06)
  call Work % Disconnect_Real_Cell(var_ins, v2_calc, kin_vis_t, phi_save)

  call Profiler % Stop('Save_Vtu_Results')

  end subroutine

!==============================================================================!
  Character(SL) function Var_Name(str1, str2, bool)
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  character(*) :: str1, str2
  logical :: bool
!==============================================================================!
  Var_Name = trim(str1)
  if (bool) Var_Name = trim(Var_Name) // " " // trim(str2)

  end function
