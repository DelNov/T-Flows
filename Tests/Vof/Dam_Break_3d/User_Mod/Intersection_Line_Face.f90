!==============================================================================!
  subroutine Intersection_Line_Face(Flow, Vof, s, pab, pint, inters)
!------------------------------------------------------------------------------!
!   This function finds intersection between a line an face s                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Vof_Type),   target :: Vof
  integer                  :: s
  real                     :: pab(2,3)
  real                     :: pint(3)
  logical                  :: inters
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, j, n
  integer,     allocatable :: list_nodes(:)
  integer,     allocatable :: edges_count(:)
  real                     :: determ, dotprod
  real                     :: lab(3), p01(3), p02(3)
  real                     :: t, u, v
  real                     :: n_unit(3), v1(3), v2(3)
  real                     :: epsloc
  real                     :: xmin, xmax, ymin, ymax
!==============================================================================!

  ! Take aliases
  Grid  => Flow % pnt_grid

  epsloc = epsilon(epsloc)

  allocate(list_nodes(Grid % faces_n_nodes(s)),   &
           edges_count(Grid % faces_n_nodes(s)))

  edges_count = 0

  do n = 1, Grid % faces_n_nodes(s)
    list_nodes(n) = Grid % faces_n(n,s)
  end do

  lab = pab(2,:) - pab(1,:)

  p01(1) = Grid % xn(list_nodes(1)) - Grid % xf(s)
  p01(2) = Grid % yn(list_nodes(1)) - Grid % yf(s)
  p01(3) = Grid % zn(list_nodes(1)) - Grid % zf(s)

  p02(1) = Grid % xn(list_nodes(2)) - Grid % xf(s)
  p02(2) = Grid % yn(list_nodes(2)) - Grid % yf(s)
  p02(3) = Grid % zn(list_nodes(2)) - Grid % zf(s)

  determ = dot_product(-lab, Math % Cross_Product(p01, p02))

  inters = .false.
  if (abs(determ) > epsloc) then
    t = dot_product( Math % Cross_Product(p01, p02), pab(1,:)       &
                 - (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/) )   &
                 / determ
    u = dot_product( Math % Cross_Product(p02, -lab), pab(1,:)      &
                 - (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/) )   &
                 / determ
    v = dot_product( Math % Cross_Product(-lab, p01), pab(1,:)      &
                 - (/Grid % xf(s), Grid % yf(s), Grid % zf(s)/) )   &
                 / determ

    if (t >= 0.0 .and. t <= 1.0) then
      pint = pab(1,:) + lab * t

      v1 = (/Grid % xn(list_nodes(2)) - Grid % xn(list_nodes(1)),   &
             Grid % yn(list_nodes(2)) - Grid % yn(list_nodes(1)),   &
             Grid % zn(list_nodes(2)) - Grid % zn(list_nodes(1))/)

      v2 = (/Grid % xn(list_nodes(3)) - Grid % xn(list_nodes(1)),   &
             Grid % yn(list_nodes(3)) - Grid % yn(list_nodes(1)),   &
             Grid % zn(list_nodes(3)) - Grid % zn(list_nodes(1))/)

      n_unit = Math % Cross_Product(v1, v2)    &
             / norm2(Math % Cross_Product(v1, v2))

      ! Verify if pint is inside face
      do n = 1, size(list_nodes)
        if (n == size(list_nodes)) then
          v1 = (/Grid % xn(list_nodes(1)) - Grid % xn(list_nodes(n)),       &
                 Grid % yn(list_nodes(1)) - Grid % yn(list_nodes(n)),       &
                 Grid % zn(list_nodes(1)) - Grid % zn(list_nodes(n))/)
        else
          v1 = (/Grid % xn(list_nodes(n + 1)) - Grid % xn(list_nodes(n)),   &
                 Grid % yn(list_nodes(n + 1)) - Grid % yn(list_nodes(n)),   &
                 Grid % zn(list_nodes(n + 1)) - Grid % zn(list_nodes(n))/)
        end if
        v2 = (/pint(1) - Grid % xn(list_nodes(n)),   &
               pint(2) - Grid % yn(list_nodes(n)),   &
               pint(3) - Grid % zn(list_nodes(n))/)
        if (dot_product(n_unit, Math % Cross_Product(v1, v2)) > 0.0) then
          edges_count(n) = 1
        end if
      end do

      if (sum(edges_count) == size(list_nodes)) then
        inters = .true.
      end if
    end if

  end if

  end subroutine
