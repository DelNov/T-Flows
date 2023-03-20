!==============================================================================!
  function Func_2(Pol, x, y, z)
!------------------------------------------------------------------------------!
!   Torus with major radius 0.2, minor radius 0.1 and centered at (0.5,0.5,0.5)!
!   Exact volume encolsed by the interface=2*PI**2*(2/3)*(1/3)**2              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  real                   :: Func_2
  real                   :: x, y, z
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Pol)
!==============================================================================!

  Func_2=(0.1)**2-((0.2)-((x-0.5)**2+(z-0.5)**2)**0.5)**2-  &
         (y-0.5)**2

  end
