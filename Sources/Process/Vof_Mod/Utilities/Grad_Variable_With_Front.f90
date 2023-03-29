!==============================================================================!
  subroutine Grad_Variable_With_Front(Vof, var, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field Flow                          !
!                                                                              !
!   (Closely related (derived from) to Field_Mod_Grad_Variable)                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type)   :: Vof
  type(Var_Type)    :: var
  real, intent(in)  :: phif
!------------------------------[Local parameters]------------------------------!
  logical, parameter ::       PRESC_GRAD = .false.
  logical, parameter :: NEW_GRAD_GIORGIA = .false.
  logical, parameter ::   NEW_GRAD_BOJAN = .false.
  integer            :: NGSZ = 48  ! number of checked neighbours
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e, c1, c2, s, c_fr, cn_l, cn_v, m
  real                      :: dist_xz, dist_xyz, x, y, z
  real                      :: alpha_xz, alpha_xyz, psi
  real                      :: x_bub, y_bub, z_bub
  real                      :: x_tr, y_tr, z_tr
  real, pointer, contiguous :: x_fr(:), y_fr(:), z_fr(:)
  integer                   :: n, cn, i_cel, i_nod, i_fac, b_fac
  real                      :: xe, ye, ze, xcn, ycn, zcn, xn, yn, zn, n_v_norm
  real, dimension(3)        :: n_v, d1, d1_n, d2, d2_n, d3, d3_n
  real, pointer, contiguous :: x_grad(:), y_grad(:), z_grad(:)
  integer                       :: j_cel
  real, dimension(3)            :: b_v
  real, pointer, contiguous     :: norm_grad_1(:), norm_grad_0(:)
  real                          :: temp_b, area_el
  real                          :: xb, yb, zb, dx, dy, dz, c_size, dg
  integer, pointer, contiguous  :: cell_list(:)
  real, pointer, contiguous     :: cell_dist(:)
!==============================================================================!

  ! Take alias
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! Calculate gradient matrix with front (reduced stencil close to front)
  call Vof % Calculate_Grad_Matrix_With_Front()

!   ! Refresh buffers for variable
!   call Grid % Exchange_Cells_Real(var % n)
!
!   ! Compute individual gradients without refreshing buffers
!   call Vof % Grad_Component_No_Refresh_With_Front(var % n, 1, var % x, phif)
!   call Vof % Grad_Component_No_Refresh_With_Front(var % n, 2, var % y, phif)
!   call Vof % Grad_Component_No_Refresh_With_Front(var % n, 3, var % z, phif)
!
!   ! Refresh buffers for gradient components
!   call Grid % Exchange_Cells_Real(var % x)
!   call Grid % Exchange_Cells_Real(var % y)
!   call Grid % Exchange_Cells_Real(var % z)

  !----------------------------------------------------------------------!
  !   PRESCRIBED TEMPERATURE GRADIENT TO TEST THE ISSUES CONCERNING      !
  !   CALCULATED GRADIENTS IN CELLS WITH FRONT. EVENTUALLY, THE SOURCE   !
  !   Vof % Grad_Component_No_Refresh_With_Front IS THE PROBLEMATIC ONE. !
  !   THE IMPOSED GRADIENT ALSO WORKS WITH RISING BUBBLE CASE (GROWTH)   !
  !----------------------------------------------------------------------!
  if(PRESC_GRAD) then
    call Work % Connect_Real_Cell(x_fr, y_fr, z_fr)
    c_fr = 0.0
    x_fr(:) = 0.0
    y_fr(:) = 0.0
    z_fr(:) = 0.0
    x_bub = 0.0
    y_bub = 0.0
    z_bub = 0.0
    
    ! Find cells within the bubble
    do c = 1, Grid % n_cells
      if(Vof % Front % elem_in_cell(c) .ne. 0) then
        c_fr = c_fr + 1
        x_fr(c_fr) =  Grid % xc(c)
        y_fr(c_fr) =  Grid % yc(c)
        z_fr(c_fr) =  Grid % zc(c)
      end if
    end do
    
    ! Approximate bubble center
    x_bub = sum(x_fr) / c_fr
    y_bub = sum(y_fr) / c_fr
    z_bub = sum(z_fr) / c_fr
    print*, "XB:", x_bub, "ZB:", z_bub
    print*, "NUMBER OF CELLS IN FRONT:", c_fr
    do c = 1, Grid % n_cells
      x = Grid % xc(c)
      y = Grid % yc(c)
      z = Grid % zc(c)
      x_tr = x + x_bub
      y_tr = y + y_bub
      z_tr = z + z_bub
   
      dist_xz = sqrt(x_tr**2 + z_tr**2)
      dist_xyz = sqrt(x_tr**2 + y_tr**2 + z_tr**2)
      alpha_xz = atan2(z_tr,x_tr)
      alpha_xyz = atan2(y_tr,x_tr)
      psi = acos(z_tr/dist_xyz)
!     if(Vof % Front % elem_in_cell(c) .eq. 0) then ! check if c is at interface
        if(Vof % fun % n(c) < 0.5) then
          ! vapor phase
          ! 2D case
          Flow % t % x(c) = 0.0
          Flow % t % z(c) = 0.0
          ! 3D case
!           Flow % t % x(c) = 0.0
!           Flow % t % y(c) = 0.0
!           Flow % t % z(c) = 0.0
        end if
        if(Vof % fun % n(c) > 0.5) then
          ! liquid phase
          ! 2D case (x-z plane)
          Flow % t % x(c) = 10000.0 * cos(alpha_xz)
          Flow % t % z(c) = 10000.0 * sin(alpha_xz)
          ! 3D case x, y, z spherical bubble
!            Flow % t % x(c) = 10000.0 * sin(psi)*cos(alpha_xyz)
!            Flow % t % y(c) = 10000.0 * sin(psi)*sin(alpha_xyz)
!            Flow % t % z(c) = 10000.0 * cos(psi)
        end if
!       end if
    end do
 
    print*, "NEW X", x_tr, "NEW Z", z_tr
    call Work % Disconnect_Real_Cell(x_fr, y_fr, z_fr)
  end if

  !----------------------------------------------------------------------!
  !   NEW IMPLEMENTATION BY GIORGIA OF THE GRADIENT INSPIRED             !
  !   BY THE WORK FROM PEREZ-RAYA ET AL. (2018). THE ALGORITHM           !
  !   WAS MODIFIED TO COMPLY TO WHAT WAS ALREADY IMPLEMENTED IN T-FLOWS. !
  !----------------------------------------------------------------------!

  if(NEW_GRAD_GIORGIA) then
    call Work % Connect_Real_Cell(x_grad, y_grad, z_grad)

    ! Browse through cells
    do c = 1, Grid % n_cells
    ! Interface element at cell c
      e = Vof % Front % elem_in_cell(c)

      ! Normal vector to element surface
      n_v = (/Vof % nx(c), Vof % ny(c), Vof % nz(c)/)
      ! Norm of normal vector
      n_v_norm = norm2(n_v)      
      cn_l = 0
      cn_v = 0
      x_grad(:) = 0.0
      y_grad(:) = 0.0
      z_grad(:) = 0.0
     
    ! There is an interface element e in cell c
      if(e > 0) then

        ! Calculate coordinates of cell center of c
        x =  Grid % xc(c)
        y =  Grid % yc(c)
        z =  Grid % zc(c)
        ! Calculate coordinates of cell center of element e
        xe = Vof % Front % Elem(e) % xe
        ye = Vof % Front % Elem(e) % ye
        ze = Vof % Front % Elem(e) % ze

        ! Define vector d1
        d1 = (/x-xe, y-ye, z-ze/)
        d1_n = dot_product(d1, n_v/n_v_norm)*n_v/n_v_norm
       
        do i_nod = 1, abs(Grid % cells_n_nodes(c)) ! local for cell
          n = Grid % cells_n(i_nod, c)             ! global
         ! Take node coordinates
          xn = Grid % xn(n)
          yn = Grid % yn(n)
          zn = Grid % zn(n)
          d2 = (/xn-xe, yn-ye, zn-ze/)
          d2_n = dot_product(d2, n_v/n_v_norm)*n_v/n_v_norm
         
          do i_cel = 1, Grid % nodes_n_cells(n)    ! local
            cn = Grid % nodes_c(i_cel, n)          ! global
           
           ! Liquid phase
            if( Vof % fun % n(c) > 0.5) then
              if(cn > 0                                .and.  &
                 cn .ne. c                             .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0 .and.  &
                 Vof % fun % n(cn) > 0.99) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)
                d3 = (/xcn-xe, ycn-ye, zcn-ze/)
                d3_n = dot_product(d3, n_v/n_v_norm)*n_v/n_v_norm
                if(norm2(d3_n) > norm2(d2_n)) then
                cn_l = cn_l + 1
                x_grad(cn) = Flow % t % x(cn)
!                y_grad(cn) = (Flow % t % n(cn) - Vof % t_sat) / d3_n(2)
                z_grad(cn) = Flow % t % z(cn)    
                 end if    
              end if
              Flow % t % x(c) = sum(x_grad) / cn_l
!              Flow % t % y(c) = sum(y_grad) / cn_l
              Flow % t % z(c) = sum(z_grad) / cn_l
            end if  
           
            if(Vof % fun % n(c) < 0.5) then
              if(cn > 0                                .and.  &
                 cn .ne. c                             .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0 .and.  &
                 Vof % fun % n(cn) < 0.5) then
                xcn = Grid % xc(cn)
                ycn = Grid % yc(cn)
                zcn = Grid % zc(cn)
                d3 = (/xcn-xe, ycn-ye, zcn-ze/)
                d3_n = dot_product(d3, n_v/n_v_norm)*n_v/n_v_norm
                if(norm2(d3_n) > norm2(d2_n)) then
                  cn_v = cn_v + 1
                  x_grad(cn) = Flow % t % x(cn)
!                  y_grad(cn) = (Flow % t % n(cn) - Vof % t_sat) / d3_n(2)
                  z_grad(cn) = Flow % t % z(cn)    
                end if      
              end if
              Flow % t % x(c) = sum(x_grad) / cn_v
!              Flow % t % y(c) = sum(y_grad) / cn_v
              Flow % t % z(c) = sum(z_grad) / cn_v
            end if
           
          end do
        end do

      end if
    end do
    call Work % Disconnect_Real_Cell(x_grad, y_grad, z_grad)
  end if
 
  !----------------------------------------------------------------------!
  !   NEW IMPLEMENTATION BY BOJAN OF THE GRADIENT INSPIRED               !
  !   BY THE WORK FROM PEREZ-RAYA ET AL. (2018).                         !
  !----------------------------------------------------------------------!

  if(NEW_GRAD_BOJAN) then
   call Work % Connect_Real_Cell(norm_grad_1, norm_grad_0, cell_dist)
   call Work % Connect_Int_Cell(cell_list)
   
   norm_grad_0(:) = 0.0
   norm_grad_1(:) = 0.0
   
   do c = 1, Grid % n_cells
      e = Vof % Front % elem_in_cell(c)
      if(e > 0) then
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
        dg = 1.5 * c_size
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
                 Vof % fun % n(cn) > 0.99 .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0) then
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
                 Vof % fun % n(cn) < 0.5 .and.  &
                 Vof % Front % elem_in_cell(cn) .eq. 0) then
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
          norm_grad_0(c) = (Vof % t_sat - temp_b) / dg
          write(50,*) xb, yb, zb, temp_b
!          Flow % t % x(c) = norm_grad_0(c) * n_v(1)
!          Flow % t % y(c) = norm_grad_0(c) * n_v(2)
!          Flow % t % z(c) = norm_grad_0(c) * n_v(3)
        end if
      end if
    end do

   call Work % Disconnect_Real_Cell(norm_grad_1, norm_grad_0, cell_dist)
   call Work % Disconnect_Int_Cell(cell_list)
  end if
 
  ! Recover the original gradient matrix
  call Flow % Calculate_Grad_Matrix()

  end subroutine
