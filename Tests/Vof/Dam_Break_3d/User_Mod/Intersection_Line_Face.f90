!==============================================================================!
  subroutine Intersection_Line_Face(flow, Vof, s, pab, pint, inters)
!------------------------------------------------------------------------------!
!   This function finds intersection between a line an face s                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Vof_Type),   target :: Vof
  integer                  :: s
  real                     :: pab(2,3)
  real                     :: pint(3)
  logical                  :: inters
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
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
  grid  => flow % pnt_grid

  epsloc = epsilon(epsloc)

  allocate(list_nodes(grid % faces_n_nodes(s)),   &
           edges_count(grid % faces_n_nodes(s)))

  edges_count = 0

  do n = 1, grid % faces_n_nodes(s)
    list_nodes(n) = grid % faces_n(n,s)
  end do

  lab = pab(2,:) - pab(1,:)

  p01(1) = grid % xn(list_nodes(1)) - grid % xf(s)
  p01(2) = grid % yn(list_nodes(1)) - grid % yf(s)
  p01(3) = grid % zn(list_nodes(1)) - grid % zf(s)

  p02(1) = grid % xn(list_nodes(2)) - grid % xf(s)
  p02(2) = grid % yn(list_nodes(2)) - grid % yf(s)
  p02(3) = grid % zn(list_nodes(2)) - grid % zf(s)

  determ = dot_product(-lab, Math_Mod_Cross_Product(p01, p02))

  inters = .false.
  if (abs(determ) > epsloc) then
    t = dot_product( Math_Mod_Cross_Product(p01, p02), pab(1,:)       &
                   - (/grid % xf(s), grid % yf(s), grid % zf(s)/) )   &
                   / determ
    u = dot_product( Math_Mod_Cross_Product(p02, -lab), pab(1,:)      &
                   - (/grid % xf(s), grid % yf(s), grid % zf(s)/) )   &
                   / determ
    v = dot_product( Math_Mod_Cross_Product(-lab, p01), pab(1,:)      &
                   - (/grid % xf(s), grid % yf(s), grid % zf(s)/) )   &
                   / determ

    if (t >= 0.0 .and. t <= 1.0) then
      pint = pab(1,:) + lab * t

      v1 = (/grid % xn(list_nodes(2)) - grid % xn(list_nodes(1)),   &
             grid % yn(list_nodes(2)) - grid % yn(list_nodes(1)),   &
             grid % zn(list_nodes(2)) - grid % zn(list_nodes(1))/)

      v2 = (/grid % xn(list_nodes(3)) - grid % xn(list_nodes(1)),   &
             grid % yn(list_nodes(3)) - grid % yn(list_nodes(1)),   &
             grid % zn(list_nodes(3)) - grid % zn(list_nodes(1))/)

      n_unit = Math_Mod_Cross_Product(v1, v2)    &
             / norm2(Math_Mod_Cross_Product(v1, v2))

      ! Verify if pint is inside face
      do n = 1, size(list_nodes)
        if (n == size(list_nodes)) then
          v1 = (/grid % xn(list_nodes(1)) - grid % xn(list_nodes(n)),       &
                 grid % yn(list_nodes(1)) - grid % yn(list_nodes(n)),       &
                 grid % zn(list_nodes(1)) - grid % zn(list_nodes(n))/)
        else
          v1 = (/grid % xn(list_nodes(n + 1)) - grid % xn(list_nodes(n)),   &
                 grid % yn(list_nodes(n + 1)) - grid % yn(list_nodes(n)),   &
                 grid % zn(list_nodes(n + 1)) - grid % zn(list_nodes(n))/)
        end if
        v2 = (/pint(1) - grid % xn(list_nodes(n)),   &
               pint(2) - grid % yn(list_nodes(n)),   &
               pint(3) - grid % zn(list_nodes(n))/)
        if (dot_product(n_unit, Math_Mod_Cross_Product(v1, v2)) > 0.0) then
          edges_count(n) = 1
        end if
      end do

      if (sum(edges_count) == size(list_nodes)) then
        inters = .true.
      end if
    end if

  end if

  end subroutine
