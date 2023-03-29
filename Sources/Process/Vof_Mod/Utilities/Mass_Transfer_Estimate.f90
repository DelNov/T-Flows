!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer, intent(in)     :: curr_dt
  integer, intent(in)     :: ini
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
  logical, parameter :: EXTRAPOLATION = .true.
  logical, parameter :: NORM_GRAD = .false.
  integer            :: NGSZ = 48  ! number of checked neighbours
  real,    parameter :: L = 2.26e+6  ! Latent heat [J/kg]   
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer     :: Grid
  type(Field_Type), pointer     :: Flow
  type(Var_Type),   pointer     :: t
  integer                       :: e, s, c1, c2, i_ele
  real                          :: cond_1, cond_2, tmp
  integer                       :: c
  character(SL)                 :: fname
  integer                       :: n, cn, i_cel, i_nod, i_fac, j_cel
  real                          :: x, y, z, xe, ye, ze, xcn, ycn, zcn
  real, dimension(3)            :: n_v, b_v
  real, pointer, contiguous     :: norm_grad_1(:), norm_grad_0(:)
  real                          :: temp_b, area_el
  real                          :: xb, yb, zb, dx, dy, dz, c_size, dg
  integer, pointer, contiguous  :: cell_list(:)
  real, pointer, contiguous     :: cell_dist(:)
!==============================================================================!
 
  !------------------------------------------------!
  !   MASS TRANSFER OLD PART WITH EXTRAPOLATION    !
  !------------------------------------------------!
  
  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow
  t    => Flow % t

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  !------------------------------------------------!
  !  Compute gradients of temperature, imposing   !
  !   saturation temperature at the interface     !
  !------------------------------------------------!
  call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  
  if(DEBUG) then
    call Grid % Save_Debug_Vtu("grad-t",                               &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
  end if

  !----------------------------------------!
  !  Compute heat flux at the interface   !
  !----------------------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(any(Vof % Front % elems_at_face(1:2,s) .ne. 0)) then

      Vof % q_int(1,s) = 0.0
      Vof % q_int(2,s) = 0.0

      do i_ele = 1, 2
        e = Vof % Front % elems_at_face(i_ele, s)
        if(e .ne. 0) then

!           Take conductivities from each side of the interface
          cond_1 = Vof % phase_cond(1)
          cond_2 = Vof % phase_cond(0)
          if(Vof % fun % n(c1) < 0.5) cond_1 = Vof % phase_cond(0)
          if(Vof % fun % n(c2) > 0.5) cond_2 = Vof % phase_cond(1)

!           Compute heat fluxes from each side of the interface
!           Keep in mind that in the Stefan's case, dot product of
!           element's surface and face's surface is positive.  If
!           VOF's definition changes, I assume this would have to
!           be adjusted accordingly.
!           Units: W/(mK) * K/m * m^2 = W
          if(  Vof % Front % elem(e) % sx * Grid % sx(s)  &                      
             + Vof % Front % elem(e) % sy * Grid % sy(s)  &                      
             + Vof % Front % elem(e) % sz * Grid % sz(s) > 0.0 ) then            
                                                                                 
               Vof % q_int(1,s) = Vof % q_int(1,s)                        &            
                    + cond_1 * ( t % x(c1) * Vof % Front % elem(e) % sx   &        
                               + t % y(c1) * Vof % Front % elem(e) % sy   &        
                               + t % z(c1) * Vof % Front % elem(e) % sz)            
                                                                                 
               Vof % q_int(2,s) = Vof % q_int(2,s)                        &        
                    + cond_2 * ( t % x(c2) * Vof % Front % elem(e) % sx   &        
                               + t % y(c2) * Vof % Front % elem(e) % sy   &        
                               + t % z(c2) * Vof % Front % elem(e) % sz)      
          else                                                                  
                                                                                 
               Vof % q_int(1,s) = Vof % q_int(1,s)                        &        
                    - cond_1 * ( t % x(c1) * Vof % Front % elem(e) % sx   &        
                               + t % y(c1) * Vof % Front % elem(e) % sy   &        
                               + t % z(c1) * Vof % Front % elem(e) % sz)            
                                                                                 
               Vof % q_int(2,s) = Vof % q_int(2,s)                        &        
                    - cond_2 * ( t % x(c2) * Vof % Front % elem(e) % sx   &        
                               + t % y(c2) * Vof % Front % elem(e) % sy   &        
                               + t % z(c2) * Vof % Front % elem(e) % sz)       
        end if


        end if  ! e .ne. 0
      end do    ! i_ele

      ! Accumulate sources in the cells surroundig the face
      ! Unit: W * kg/J = kg / s
!      tmp = sqrt((Vof % Front % Elem(e) % sx)**2    &
!               + (Vof % Front % Elem(e) % sy)**2    &
!               + (Vof % Front % Elem(e) % sz)**2)
!               
!      if(Vof % Front % elem_in_cell(c1) .ne. 0) then
!        Vof % m_dot(c1) =  Vof % m_dot(c1)                       &
!                        + (Vof % q_int(2,s) - Vof % q_int(1,s))  &
!                        / 2.26e+6
!        write(100,*) c1, Vof % m_dot(c1), tmp
!      end if
!
!      if(Vof % Front % elem_in_cell(c2) .ne. 0) then
!        Vof % m_dot(c2) =  Vof % m_dot(c2)                       &
!                        + (Vof % q_int(2,s) - Vof % q_int(1,s))  &
!                        / 2.26e+6
!        write(101,*) c2, Vof % m_dot(c2), tmp
!      end if

    end if      ! face is at the front

  end do

  !-------------------------------------------------!
  !   MASS TRANSFER WITH GRADIENT EXTRAPOLATION     !
  !-------------------------------------------------!
  if(EXTRAPOLATION) then
  ! Gradients are fresh here
  ! Intialize t_0 and t_1
    Vof % t_0 % x(1:) = Flow % t % x(1:)
    Vof % t_0 % y(1:) = Flow % t % y(1:)
    Vof % t_0 % z(1:) = Flow % t % z(1:)
    Vof % t_1 % x(1:) = Flow % t % x(1:)
    Vof % t_1 % y(1:) = Flow % t % y(1:)
    Vof % t_1 % z(1:) = Flow % t % z(1:)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % x, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % y, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % z, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % x, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % y, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % z, towards=0)
  
    if(DEBUG) then
      fname = "grad-t-0-"
      write(fname(10:15), '(i5.5)') curr_dt
      call Grid % Save_Debug_Vtu(fname,                            &
                                 scalar_cell = Vof % fun % n,      &
                                 scalar_name = "vof_fun",          &
                                 vector_cell = (/Vof % t_0 % x,    &
                                                 Vof % t_0 % y,    &
                                                 Vof % t_0 % z/),  &
                                 vector_name = "grad-t-0")
      fname = "grad-t-1-"
      write(fname(10:15), '(i5.5)') curr_dt
      call Grid % Save_Debug_Vtu(fname,                            &
                                 scalar_cell = Vof % fun % n,      &
                                 scalar_name = "vof_fun",          &
                                 vector_cell = (/Vof % t_1 % x,    &
                                                 Vof % t_1 % y,    &
                                                 Vof % t_1 % z/),  &
                                 vector_name = "grad-t-1")
    end if
  
    Vof % q_0(:) = 0
    Vof % q_1(:) = 0
    do c = 1, Grid % n_cells
      e = Vof % Front % elem_in_cell(c)
      if(e > 0) then
        ! Heat flux to the interface in cell c from phase 0
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_0(c) = (  Vof % t_0 % x(c) * Vof % Front % Elem(e) % sx   &
                        + Vof % t_0 % y(c) * Vof % Front % Elem(e) % sy   &
                        + Vof % t_0 % z(c) * Vof % Front % Elem(e) % sz)  &
                     * Vof % phase_cond(0)

        ! Heat flux to the interface in cell c from phase 1
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_1(c) = (  Vof % t_1 % x(c) * Vof % Front % Elem(e) % sx   &
                        + Vof % t_1 % y(c) * Vof % Front % Elem(e) % sy   &
                        + Vof % t_1 % z(c) * Vof % Front % Elem(e) % sz)  &
                     * Vof % phase_cond(1)
        ! Units: W / (J/kg) = kg/s
        Vof % m_dot(c) = (Vof % q_1(c) - Vof % q_0(c)) / L
        
        tmp = sqrt((Vof % Front % Elem(e) % sx)**2 + (Vof % Front % Elem(e) % sy)**2 + (Vof % Front % Elem(e) % sz)**2)
        write(100,*) c, Vof % m_dot(c), tmp
      end if
    end do
  end if 
  
  !------------------------------------------------------!
  !   MASS TRANSFER WITH NORMAL GRADIENT CALCULATION     !
  !------------------------------------------------------!
  if(NORM_GRAD) then
  
   call Work % Connect_Real_Cell(norm_grad_1, norm_grad_0, cell_dist)
   call Work % Connect_Int_Cell(cell_list)
   
   norm_grad_0(:) = 0.0
   norm_grad_1(:) = 0.0
   
   do c = 1, Grid % n_cells
      e = Vof % Front % elem_in_cell(c)
!      if(e > 0) then
        x =  Grid % xc(c)
        y =  Grid % yc(c)
        z =  Grid % zc(c)
        xe = Vof % Front % Elem(e) % xe
        ye = Vof % Front % Elem(e) % ye
        ze = Vof % Front % Elem(e) % ze

        ! Unity normal vector
        n_v = (/Vof % nx(c), Vof % ny(c), Vof % nz(c)/)

        ! Work out the cell size locally
        c_size = HUGE
        do i_fac = 1, Grid % cells_n_faces(c)
          dx = abs(Grid % dx(i_fac))
          dy = abs(Grid % dy(i_fac))
          dz = abs(Grid % dz(i_fac))
          c_size = min(c_size, sqrt((dx**2 + dy**2 + dz**2)))
        end do
        dg = 2.0 * c_size
        b_v = n_v * dg
        xb = xe + b_v(1)
        yb = ye + b_v(2)
        zb = ze + b_v(3)
       
       ! Liquid phase    
        if(Vof % fun % n(c) > 0.5) then
          cell_dist(:) = HUGE
          cell_list(:) = 0
          ! Calculate the distances of all cells around
          ! the node and store them in cell_dist (and cell_list)
          j_cel = 0
          do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local
            n = Grid % cells_n(i_nod, c)             ! global
            do i_cel = 1, Grid % nodes_n_cells(n)    ! local
              cn = Grid % nodes_c(i_cel, n)          ! global

              if(cn > 0                   .and.  &
                 cn .ne. c                .and.  &
                 Vof % fun % n(cn) > 0.99) then
!                 Vof % Front % elem_in_cell(cn) .eq. 0) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)
                
                j_cel = j_cel + 1

                cell_dist(j_cel) = sqrt( (xb - xcn)**2     &
                                       + (yb - ycn)**2     &
                                       + (zb - zcn)**2 )
                cell_list(j_cel) = cn
              end if
            end do
          end do

          ! Find out the cell with shortest distance to b
          call Sort % Real_Carry_Int(cell_dist(1:j_cel), cell_list(1:j_cel))

          ! Closest cell's index is now in cell_list(1)
          cn  = cell_list(1)
          xcn = Grid % xc(cn)
          ycn = Grid % yc(cn)
          zcn = Grid % zc(cn)
          temp_b = Flow % t % n(cn) + (Flow % t % x(cn) * (xb - xcn) &
                                     + Flow % t % y(cn) * (yb - ycn) &
                                     + Flow % t % z(cn) * (zb - zcn))
                                     
          ! Calculate normal gradient in liquid phase
          norm_grad_1(c) = (temp_b - Vof % t_sat) / dg
          write(40,*) xb, yb, zb, temp_b
!          Flow % t % x(c) = norm_grad_1(c) * n_v(1)
!          Flow % t % y(c) = norm_grad_1(c) * n_v(2)
!          Flow % t % z(c) = norm_grad_1(c) * n_v(3)
        end if
        
       ! Vapor phase
        if(Vof % fun % n(c) < 0.5) then
          n_v = - n_v
          b_v = n_v * dg
          xb = xe + b_v(1)
          yb = ye + b_v(2)
          zb = ze + b_v(3)
          cell_dist(:) = HUGE
          cell_list(:) = 0

          ! Calculate the distances of all cells around
          ! the node and store them in cell_dist (and cell_list)
          j_cel = 0
          do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local
            n = Grid % cells_n(i_nod, c)             ! global
            do i_cel = 1, Grid % nodes_n_cells(n)    ! local
              cn = Grid % nodes_c(i_cel, n)          ! global

              if(cn .ne. c               .and.  &
                 Vof % fun % n(cn) < 0.5) then
!                 Vof % Front % elem_in_cell(cn) .eq. 0) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)
                
                j_cel = j_cel + 1
                
                cell_dist(j_cel) = sqrt(  (xb - xcn)**2  &
                                        + (yb - ycn)**2  &
                                        + (zb - zcn)**2 )
                cell_list(j_cel) = cn
              end if
            end do
          end do

          ! Find out the cell with shortest distance to b
          call Sort % Real_Carry_Int(cell_dist(1:j_cel), cell_list(1:j_cel))

          ! Closest cell's index is now in cell_list(1)
          cn  = cell_list(1)
          xcn = Grid % xc(cn)
          ycn = Grid % yc(cn)
          zcn = Grid % zc(cn)
          temp_b = Flow % t % n(cn) + (Flow % t % x(cn) * (xb - xcn) &
                                     + Flow % t % y(cn) * (yb - ycn) &
                                     + Flow % t % z(cn) * (zb - zcn))
          ! Calculate normal gradient in vapor phase
          norm_grad_0(c) = (temp_b - Vof % t_sat) / dg
          write(50,*) xb, yb, zb, temp_b
!          Flow % t % x(c) = norm_grad_0(c) * n_v(1)
!          Flow % t % y(c) = norm_grad_0(c) * n_v(2)
!          Flow % t % z(c) = norm_grad_0(c) * n_v(3)
        end if
!      end if
    end do
    
    Vof % q_0(:) = 0.0
    Vof % q_1(:) = 0.0
    do c = 1, Grid % n_cells
      e = Vof % Front % elem_in_cell(c)
      if(e > 0) then
        ! Total front element area
        area_el = sqrt((Vof % Front % Elem(e) % sx)**2    &
                     + (Vof % Front % Elem(e) % sy)**2    &
                     + (Vof % Front % Elem(e) % sz)**2)
        ! Heat flux to the interface in cell c from phase 0
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_0(c) = Vof % phase_cond(0) * norm_grad_0(c) * area_el
!        if(Vof % fun % n(c) < 0.5) then
!          Vof % q_0(c) = (  Flow % t % x(c) * Vof % Front % Elem(e) % sx   &
!                          + Flow % t % y(c) * Vof % Front % Elem(e) % sy   &
!                          + Flow % t % z(c) * Vof % Front % Elem(e) % sz)  &
!                       * Vof % phase_cond(0)
!        end if
        ! Heat flux to the interface in cell c from phase 1
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_1(c) = Vof % phase_cond(1) * norm_grad_1(c) * area_el
!        if(Vof % fun % n(c) > 0.5) then
!          Vof % q_1(c) = (  Flow % t % x(c) * Vof % Front % Elem(e) % sx   &
!                          + Flow % t % y(c) * Vof % Front % Elem(e) % sy   &
!                          + Flow % t % z(c) * Vof % Front % Elem(e) % sz)  &
!                       * Vof % phase_cond(1)
!        end if

        ! Units: W / (J/kg) = kg/s
        Vof % m_dot(c) = (Vof % q_1(c) - Vof % q_0(c)) / L
        write(100,*) c, Vof % m_dot(c), area_el
        write(200,*) c, norm_grad_1(c), norm_grad_0(c)
      end if
    end do
    
    if(DEBUG) then
      fname = "vap-grad-"
      write(fname(10:15), '(i5.5)') curr_dt
      call Grid % Save_Debug_Vtu(fname,                            &
                                 scalar_cell = norm_grad_0,        &
                                 scalar_name = "grad-norm-0")

      fname = "liq-grad-"
      write(fname(10:15), '(i5.5)') curr_dt
      call Grid % Save_Debug_Vtu(fname,                            &
                                 scalar_cell = norm_grad_1,        &
                                 scalar_name = "grad-norm-1")
    end if
 
   call Work % Disconnect_Real_Cell(norm_grad_1, norm_grad_0, cell_dist)
   call Work % Disconnect_Int_Cell(cell_list)
  end if 
  
  end subroutine
