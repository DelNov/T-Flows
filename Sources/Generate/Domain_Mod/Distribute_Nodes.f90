!==============================================================================!
  subroutine Distribute_Nodes(Dom, Grid, b, bw,   &
                              is, js, ks, ie, je, ke)
!------------------------------------------------------------------------------!
!>  This subroutine places nodes on a line defined by local block positions.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initialization: Retrieves the resolution of the block (ni, nj, nk) and   !
!     the starting coordinates (x0, y0, z0).                                   !
!   * Calculating delta values: Determines the change in coordinates (delx,    !
!     dely, delz) across the specified range within the block.                 !
!   * Early exit condition: If the number of nodes to distribute is less       !
!     than 2, the subroutine returns immediately, as no distribution is        !
!     needed.                                                                  !
!   * Node distribution logic:                                                 !
!     - Linear distribution (positive weight): If the weight (bw) is positive, !
!       it distributes the nodes linearly between the start and end points,    !
!       adjusting the node positions based on a computed delta shift (dt).     !
!     - Hyperbolic distribution (negative weight): For negative weights, it    !
!       uses a hyperbolic tangent function (tanh) to distribute the nodes.     !
!       This involves determining a case based on the weight value and         !
!       adjusting the node positions accordingly.                              !
!   * Iterating over block dimensions: The subroutine iterates over the range  !
!     specified by is, js, ks to ie, je, ke within the block, placing the      !
!     nodes according to the calculated distribution.                          !
!   * Updating node coordinates: For each position in the range, it updates    !
!     the node coordinates (xn, yn, zn) in the Grid based on the selected      !
!     distribution method (linear or hyperbolic).                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom         !! domain in which the grid is generated
  type(Grid_Type)     :: Grid        !! grid being generated
  integer, intent(in) :: b           !! current block
  integer, intent(in) :: is, js, ks  !! starting index in a logical direction
  integer, intent(in) :: ie, je, ke  !! ending index in a logical direction
  real,    intent(in) :: bw          !! weight of the block
!----------------------------------[Calling]-----------------------------------!
  real :: atanh
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, ni, nj, nk, i, j, k, node, case
  real    :: x0, y0, z0, delx, dely, delz, t, dt, ddt, pr, xi
!==============================================================================!

  ni = Dom % blocks(b) % resolutions(1)
  nj = Dom % blocks(b) % resolutions(2)
  nk = Dom % blocks(b) % resolutions(3)

  x0   = Grid % xn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  y0   = Grid % yn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  z0   = Grid % zn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  delx = Grid % xn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - x0
  dely = Grid % yn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - y0
  delz = Grid % zn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - z0

  n = max( (ie-is), (je-js),  (ke-ks) )

  if(n < 2) return

  !-------------------------!
  !   Linear distribution   !
  !-------------------------!
  if(bw > 0.0) then
    ddt = ( 2.0*(1.0-bw) ) / ( real(n)*(real(n)-1.0)*(1.0+bw) )
    t=0.0
    node = Grid % n_nodes + ni*nj*nk  ! estimated last node
    call Grid % Allocate_Nodes(node)  ! expand nodes
    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie .ne. is ) then
            dt=1.0/real(n)+(real(i)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i+1
            if(i < ie) then
              if(.not. Grid % node_positioned(node)) then
                Grid % xn(node) = x0 + t*delx
                Grid % yn(node) = y0 + t*dely
                Grid % zn(node) = z0 + t*delz
              end if
            end if
          end if
          if( je .ne. js ) then
            dt=1.0/real(n)+(real(j)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-1)*ni*nj + (j-0)*ni + i
            if(j < je) then
              if(.not. Grid % node_positioned(node)) then
                Grid % xn(node) = x0 + t*delx
                Grid % yn(node) = y0 + t*dely
                Grid % zn(node) = z0 + t*delz
              end if
            end if
          end if
          if( ke .ne. ks ) then
            dt=1.0/real(n)+(real(k)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-0)*ni*nj + (j-1)*ni + i
            if(k < ke) then
              if(.not. Grid % node_positioned(node)) then
                Grid % xn(node) = x0 + t*delx
                Grid % yn(node) = y0 + t*dely
                Grid % zn(node) = z0 + t*delz
              end if
            end if
          end if
        end do
      end do
    end do

  !-----------------------------!
  !   Hyperbolic distribution   !
  !-----------------------------!
  else
    case = 0
    if     ((bw  >  -0.5).and.(bw <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(bw))
      case = 1
    else if((bw >=  -0.75).and.(bw  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(bw))
      case = 2
    else
      pr = -bw
      case = 3
    end if

    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie .ne. is ) then
            if(case .eq. 1) xi =       - 1.0*real(i)/real(n)
            if(case .eq. 2) xi =   1.0 - 1.0*real(i)/real(n)
            if(case .eq. 3) xi = - 1.0 + 2.0*real(i)/real(n)
            node = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i+1
            if( (i < ie).and.(.not. Grid % node_positioned(node)) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
          if( je .ne. js ) then
            if(case .eq. 1) xi =        - 1.0*real(j)/real(n)
            if(case .eq. 2) xi =    1.0 - 1.0*real(j)/real(n)
            if(case .eq. 3) xi = -  1.0 + 2.0*real(j)/real(n)
            node = Grid % n_nodes + (k-1)*ni*nj + (j-0)*ni + i
            if( (j < je).and.(.not. Grid % node_positioned(node)) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
          if( ke .ne. ks ) then
            if(case .eq. 1) xi =       - 1.0*real(k)/real(n)
            if(case .eq. 2) xi =   1.0 - 1.0*real(k)/real(n)
            if(case .eq. 3) xi = - 1.0 + 2.0*real(k)/real(n)
            node = Grid % n_nodes + (k-0)*ni*nj + (j-1)*ni + i
            if( (k < ke).and.(.not. Grid % node_positioned(node)) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % yn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
        end do
      end do
    end do

  end if

  end subroutine
