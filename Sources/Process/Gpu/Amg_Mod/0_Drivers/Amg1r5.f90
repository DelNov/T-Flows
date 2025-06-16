!==============================================================================!
  subroutine Amg1r5(Amg, a_val, a_row, a_col, phi, b, n, eps, iswtch)
!------------------------------------------------------------------------------!
!
!   Amg1r5
!
!   Release 1.5, October 1990
!
!   Changes against version 1.1, july 1985:
!
!   1. A bug was detected which under certain circumstances influenced
!      slightly the convergence rate of Amg1r1.  For that reason, the
!      following line in subroutine resc (BN: Restrict_Residuals):
!      iw(imaxw(kc-1)+1) = ia(imin(kc)) has been changed to:
!      iw(imaxw(kc-1)+1) = iaux.  BN: This line now reads:
!      iw(Amg % imaxw(level_c-1)+1) = iaux1
!
!   2. A bug was detected in subroutine pwint (BN: Interpolation_Weights.)
!      Under certain circumstances an undefined variable was used.  Although
!      this did not affect the numerical results, problems can occur if
!      checking for undefined variables is used. To fix this error, in pwint
!      the label 1000 was moved to the statement:
!      iblck1 = Amg % iminw(level).
!
!   3. A parameter lratio has been introduced, denoting the ratio of space
!      occupied by a double precision real variable and that of an integer.
!      For the IBM-version lratio has been set to 2.  Change this value if
!      necessary. (if, for example, you want to change the double precision
!      vectors to single precision, lratio has to be set to 1. In the
!      Yale SMP - routine ndrv there is a parameter lratio, too.
!      BN: OK, this is nowadays, when we have dynamic memory allocation and
!      when we don't re-use the same memory space by integers and reals,
!      rather obsolete.  Although I have introduced a parameter to store it
!      called AMG_L_RATIO, it is a candidate for deletion, really.
!
!   4. Type declarations real*4 and real*8 have been changed to the
!      standard-conforming keywords real and double precision, respectively.
!
!   5. Calls to the following intrinsic functions have been replaced by
!      calls using generic names: dsqrt, min0, max0, iabs, dabs, float,
!      dfloat, dmax1, isign, idint, dlog10.
!
!   6. A save statement has been inserted in all subroutines.  BN: This could
!      also be one of the things of the past, particularly when Amg became
!      a Fortran module, but I keep it nonetheless as a precaution.
!
!------------------------------------------------------------------------------!
!
!   Change against version 1.3, april 1986:
!
!   1. A bug in subroutine check (BN: Check_Matrix_Properties) has been
!      removed.  If the original matrix was stored in an unsymmetric way,
!      the symmetrization by amg1r3 could fail under certain circumstances.
!      For a fix, the following statements in subroutine check have been
!      changed:
!
!      do 450 j=ia(i)+1,ia(i+1)-1 was changed to
!      do 450 j=ia(i)+1,icg(i)-1
!
!      do 430 j1=ia(i1)+1,ia(i1+1)-1 was changed to
!      do 430 j1=ia(i1)+1,icg(i1)-1
!
!      do 550 j=ia(i)+1,ia(i+1)-1 was changed to
!      do 550 j=ia(i)+1,icg(i)-1
!
!      do 530 j1=ia(i1)+1,ia(i1+1)-1 was changed to
!      do 530 j1=ia(i1)+1,icg(i1)-1
!
!  2. The explanatory part in subroutine Amg1r5 has been enlarged to avoid
!     misunderstandings in the definition of the argument list.  BN: The
!     argument list has changed so much, that this hardly matters any more.
!
!------------------------------------------------------------------------------!
!
!   Change against version 1.4, october, 1990 (by john w. ruge)
!
!   1. A bug in subroutine check (BN: Check_Matrix_Properties) has been
!      removed.  If the original matrix was stored in an unsymmetric way,
!      the symmetrization by amg1r3 could still fail under certain circum-
!      stances, and was not fixed in the previous version.  In addition, the
!      routine was changed in order to avoid some unnecessary row searches
!      for transose entries.
!
!==============================================================================!
!
!   AMG module for solving linear systems A*u=f
!
!   Assumptions on A:
!
!     The program requires:
!
!     - diagonal entries are always positive (on all grids);
!     - A is a square matrix which is either regular or singular
!       with rowsums=0.
!
!     For theoretical reasons the following should hold:
!
!     - A positive definite (or semi-definite with rowsum=0)
!     - A "essentially" positive type, i.e.,
!
!       -- diagonal entries must be > 0 ;
!       -- most of the off-diagonal entries <= 0 ;
!       -- rowsums should be >= 0 .
!
!     The user has to provide the matrix A, the right hand side f and
!     pointer vectors ia and ja.  (BN: row pointersi ia and columns ja)
!
!   --------------------------------------------------------------
!
!   Storage of A (BN: Description of CSR format):
!
!     The non-zero entries of the matrix A are stored in "compressed" sky-line
!     fashion in a 1-d vector a, i.e., row after row, each row starting with
!     its diagonal element. The other non-zero row entries follow their
!     diagonal entry in any order.  (BN: Nowadays this would be called CRS
!     format, for Compressed Row Storage.)
!
!     In order to identify each element in A, the user has to provide two
!     pointer arrays ia and ja.  If n_unknowns denotes the total number of
!     unknowns, the non-zero entries of any row i of A (1 <= i <= n_unknowns)
!     are stored in A(j) where the range of j is given by:
!
!     ia(i) .le. j .le. ia(i+1)-1
!
!     Thus, ia(i) points to the position of the diagonal entry of row i
!     within the vector A.  In particular,
!
!     ia(1) = 1 ,  ia(n_unknowns+1) = 1 + n_nonzeros
!
!     where n_nonzeros denotes the total number of matrix entries stored. The
!     pointer vector ja has to be defined such that any entry a(j) corresponds
!     to the unknown u(ja(j)), i.e., ja(j) points to the column index of a(j).
!     In particular, a(ia(i)) is the diagonal entry of row i and corresponds
!     to the unknown u(i): ja(ia(i))=i.
!
!     In this terminology, the i-th equation reads as follows
!
!     (for any i with  1.le.i.le.n_unknowns):
!
!       f(i) =    sum    a(j) * u(ja(j))
!              j1<=j<=j2
!
!     where f(i) denotes the i-th component of the right hand side and
!
!       j1 = ia(i) ,  j2 = ia(i+1)-1.
!
!     Notes:
!
!     - The entry ia(n_unknowns+1) has to point to the first free entry in
!       vectors a and ja, respectively. otherwise, Amg cannot know the length
!       of the last matrix row.
!
!     - The input vectors a, ia and ja are changed by Amg1r5.  So, after
!       return from Amg1r5, the package must not be called a second time
!       without having newly defined the input vectors and using iswtch=4.
!       Otherwise, the setup phase will fail.
!
!     - On the other hand, running Amg a second time on the same input data
!       with iswtch=4 has no sense, because the results of the first setup
!       phase are still stored and thus this phase can be skipped in a second
!       call.  In order to do this, set iswtch to 1, 2 or 3.
!
!------------------------------------------------------------------------------!
!
!   Input via arrays (see above):
!
!     a  - matrix A
!
!     ia - pointer vector (BN: row pointer)
!
!     ja - pointer vector (BN: column)
!
!     u  - first approximation to solution
!
!     f  - right hand side
!
!------------------------------------------------------------------------------!
!
!   Scalar input parameters of Amg1r5 (BN: Some of them are now module
!   member data, some are still local to Amg1r5i, some local to other
!   subroutines):
!
!   The input parameters of Amg1r5 in the list below are arranged according
!   to their importance to the general user.  The parameters preceeded by
!   a "*" must be specified explicitely.  All the other parameters are set to
!   standard values if zero on input.
!
!   There are four classes of input parameters with decreasing priority:
!
!   1. Parameters describing the user-defined problem and dimensioning of
!      vectors in the calling program
!
!   2. Parameters specifying some general algorithmic alternatives and the
!      amount of output during solution
!
!   3. Parameters controlling the multigrid cycling during the solution phase
!
!   4. Parameters controlling the creation of coarser grids and interpolation
!      formulas.
!
!   Only the Class 1 - parameters must be specified explicitely by the user.
!
!   Class 2 - parameters control the general performance of Amg1r5. Changing
!   them doesn't require understanding the Amg - algorithm.
!
!   Specifying non-standard-values for Class 3 - parameters presupposes a
!   general knowledge of multigrid methods.
!
!   Function of class 4 - parameters is only understandable after studying
!   the Amg-algorithm in detail.
!
!   Fortunately in most cases the choice of Class 3 and 4 - parameters isn't
!   critical and using the Amg1r5 - supplied standard values should give
!   satisfactory results.
!
!   --------------------------------------------------------------
!
!   Class 1 - parameters:
!
!   * nda        - Dimensioning of vector A in calling program
!
!   * ndu        - Dimensioning of vector u, f and ia in calling program
!
!   * ndw        - Dimensioning of vector kwork in calling program
!
!   BN: Since the transition to Fotran 90, nda and ndu are a bit redundant,
!   since Fortran know the dimensions of its arrays.  Thus, they are local
!   variables, where needed.  The work vector (kwork, used to be ig, I think)
!   was split into three different vectors, thus the usage of ndw will be
!   decomissioned in the future.  TLDR: these three parameters are obsolete.
!
!   * n_unknowns - Number of unknowns
!
!   * matrix     - Integer value containing info about the matrix A.
!
!     1st digit of matrix  --  isym:
!       =1: A is symmetric;
!       =2: A is not symmetric.
!
!     2nd digit of matrix  --  irow0:
!       =1: l has rowsum zero;
!       =2: l does not have rowsum zero.
!
!     BN: in order to avoig ghost numbers, I have introduced four constants
!     in module Amg, more preciselly in file Amg.h90:
!
!       AMG_SINGULAR_MATRIX      = 1,  &  ! matrix is singular (like pressure)
!       AMG_NON_SINGULAR_MATRIX  = 2,  &  ! matrix is non-singular
!       AMG_SYMMETRIC_MATRIX     = 1,  &  ! matrix is symmetric
!       AMG_NON_SYMMETRIC_MATRIX = 2      ! matrix is non-symmetric
!
!   --------------------------------------------------------------
!
!   Class 2 - parameters:
!
!     iswtch - Parameter controlling which modules of Amg1r5 are to be used.
!                =1:   call for -----, -----, -----, Wrkcnt.
!                =2:   call for -----, -----, Solve, Wrkcnt.
!                =3:   call for -----, First, Solve, Wrkcnt.
!                =4:   call for Setup, First, Solve, Wrkcnt.
!
!                Setup defines the operators needed in the solution phase.
!
!                First (BN: now First_Guess) initializes the solution vector
!                  (see parameter ifirst).
!
!                Solve computes the solution by Amg cycling
!                  (see parameter ncyc).
!
!                Wrkcnt provides the user with information about residuals,
!                  storage requirements and cp-times (see also Amg % iout).
!
!                If Amg1r5 is called the first time, iswtch has to be = 4.
!                Independent of iswtch, single modules can be bypassed by a
!                proper choice of the corresponding parameter.
!
!     Amg % iout - Parameter controlling the amount of output during solution
!              phase:
!                  =0: no output (except for messages)
!                  =1: residual before and after solution process
!                  =2: add.: statistics on cp-times and storage requirements
!                  =3: add.: residual after each Amg-cycle
!                  =4: add.: residual after each CG/BiCG solver sweep
!
!   --------------------------------------------------------------
!
!   Class 3 - parameters:
!
!     levelx - Maximum number of mg-levels to be created (>=1).
!
!     ifirst - Parameter for first approximation.
!
!                1st digit of ifirst: not used; has to be non-zero.
!
!                2nd digit of ifirst  --  itypu:
!                  =0: no setting of first approximation,
!                  =1: first approximation constant to zero,
!                  =2: first approximation constant to one,
!                  =3: first approximation is random function with the
!                      concrete random sequence being determined by the
!                      follwing digits.
!
!                rest of ifirst  --  rndu:
!                  determines the concrete random sequence used in the case
!                  itypu=3. (ifirst=13 is equivalent to ifirst=1372815)
!
!     ncyc   - Integer parameter describing the type of cycle to be
!              used and the number of cycles to be performed.
!
!                1st digit of ncyc  --  igam:
!                  =1: v -cycle,
!                  =2: v*-cycle,
!                  =3: f -cycle,
!                  =4: w -cycle.
!
!                If ncyc is negativ, then the approximation of the problem on
!                the second finest grid is computed by igam v-cycles on that
!                particular grid.
!
!                2nd digit of ncyc  --  icgr:
!                  =0: no conjugate gradient,
!                  =1: conjugate gradient (only first step of cg),
!                  =2: conjugate gradient (full cg).
!
!                3rd digit of ncyc  --  iconv: Convergence criterion for the
!                user-defined problem (finest grid):
!                  =1: perform a fixed number of cycles as given by ncycle
!                      (see below)
!                  =2: stop, if  ||res|| < eps
!                  =3: stop, if  ||res|| < eps * |f|
!                  =4: stop, if  ||res|| < eps * |u| * |diag|
!                  with ||res|| = l2-norm of residual,
!                         eps     (see input parameter eps)
!                         |f|   = supremum norm of right hand side
!                         |u|   = supremum norm of solution
!                       |diag|  = maximal diagonal entry in matrix A
!                  Note that in any case the solution process stops after at
!                  most ncycle cycles.
!
!                  rest of ncyc  --  ncycle:
!                    maximal number of cycles to be performed (>0) or
!                    ncycle=0: no cycling.
!
!     eps    - Convergence criterion for solution process: (see parameter
!              ncyc).  Note that no more than ncycle cycles are performed,
!              regardless of eps.
!
!     madapt - Integer value specifying the choice of coarsest grid in cycling:
!
!                1st digit of madapt  --  msel:
!                  =1: In cycling, all grids constructed in the setup phase
!                      are used without check.
!                  =2: The number of grids is automatically reduced if the
!                      convergence factor on the coarser grids is found to be
!                      larger than a given value fac (see below).
!
!                Rest of madapt  --  fac: The rest of madapt defines the
!                fractional part of a real number fac between 0.1 and 0.99,
!                e.g. madapt=258 means msel=2 and fac=0.58.  If madapt
!                consists of only one digit, fac is set to 0.7 by default.
!
!     def_relax_down - Parameter describing relaxation (downwards):
!
!                1st digit of def_relax_down: Not used; has to be non-zero.
!
!                2nd digit of def_relax_down  --  n_relax_down: Actual number
!                    of smoothing steps to be performed the type of which is
!                    given by the following digits
!
!                following digits  --  array type_relax_down:
!                  =1: relaxation over the f-points only
!                      (AMG_RELAX_F_POINTS)
!                  =2: full gs sweep
!                      (AMG_RELAX_FULL_GS)
!                  =3: relaxation over the c-points only
!                      (AMG_RELAX_C_POINTS)
!                  =4: full more color sweep, highest color first
!                      (AMG_RELAX_MULTICOLOR)
!
!     def_coarse_solver - Parameter controlling the solution on coarsest grid:
!
!                1st digit  --  coarse_solver:
!                  =1: Gauss-Seidel method
!                      (AMG_SOLVER_GS)
!                  =2: Conjugate Gradient solver
!                      (AMG_SOLVER_CG)
!                  =3: Bi-Conjugate Gradient solver
!                      (AMG_SOLVER_BICG)
!
!                Rest of def_coarse_solver  --  n_relax_coarse: (only if
!                coarse_solver=1) number of GS sweeps on coarsest grid (>=0).
!                If n_relax_coarse=0, then as many GS sweeps are performed as
!                are needed to reduce the residual by two orders of magnitude.
!                (maximal 100 relaxation sweeps)
!
!     def_relax_up - Parameter for relaxation (up), analogous to def_relax_down.
!
!   --------------------------------------------------------------
!
!   Class 4 - parameters:
!
!     ecg1, - Real parameters affecting the creation of coarser grids and/or
!     ecg2,   the definition of the interpolation.  The choice of these
!     ewt2    parameters depends on the actual Amg version (see subroutine
!             Coarsening)
!
!     nwt   - Integer parameter affecting the creation of coarser grids and/or
!             the definition of the interpolation.  The choice of this
!             parameter depends on the actual Amg version (see subroutine
!             Coarsening)
!
!     ntr   - Parameter controlling coarse-grid operator truncation
!               =0: pairs of zeroes are removed from coarse grid operators
!               =1: no coarse-grid operator truncation
!
!------------------------------------------------------------------------------!
!
!   Output:
!
!     u    - Contains the computed solution
!
!
!     ierr - Error parameter:
!
!              > 0: fatal error (abnormal termination of Amg1r5)
!              < 0: non-fatal error (execution of Amg1r5 con tinues)
!
!              Error codes in detail:
!
!              1. Dimensioning too small for vector
!
!                   a      (ierr = 1)
!                   ia     (ierr = 2)
!                   ja     (ierr = 3)
!                   u      (ierr = 4)
!                   f      (ierr = 5)
!                   kwork  (ierr = 6)
!                   no yale-smp because of storage (nda too small): (ierr = -1)
!                   no yale-smp because of storage (nda too small): (ierr = -3)
!                   no cg because of storage (ndu too small):       (ierr = -4)
!                   no space for transpose of interpolation (nda or nda too
!                   small): (ierr = -1)
!
!              2. Input data erroneous:
!
!                   a-entry missing, isym = 1:           (ierr = -11)
!                   parameter matrix may be erroneous:   (ierr = -12)
!                   diagonal element not stored first:   (ierr =  13)
!                   diagonal element not positiv:        (ierr =  14)
!                   pointer ia erroneous:                (ierr =  15)
!                   pointer ja erroneous:                (ierr =  16)
!                   parameter iswtch erroneous:          (ierr =  17)
!                   parameter levelx erroneous:          (ierr =  18)
!
!              3. Errors of the Amg1r5-system (should not occur):
!
!                   transpose a-entry missing:           (ierr =  21)
!                   interpolation entry missing:         (ierr =  22)
!
!              4. Algorithmic errors:
!
!                   cg-correction not defined:           (ierr =  31)
!                   no yale-smp because of error in
!                   factorization:                       (ierr = -32)
!
!------------------------------------------------------------------------------!
!
!   Work space:
!
!     The integer vector kwork has to be passed to Amg1r5 as work space.
!     BN: OK, that's how it once was, but I split the entire vector kwork
!     into three vectors: iw, ifg and icg, and pass these to subroutines
!     as work space.
!
!------------------------------------------------------------------------------!
!
!   Dimensioning of input vectors and work space:
!
!     It's impossible to tell in advance the exact storage requirements of
!     Amg.  Thus, the following formulas give only reasonable guesses for the
!     vector lengths which have to be declared in the calling program.  In
!     these formulas n_nonzeros denotes the number of non-zero entries in the
!     input-matrix A and n_unknowns is the number of unknowns.
!
!     Vector         Needed length (guess)
!       a               3   * n_nonzeros + 5 * n_unknowns
!       ja              3   * n_nonzeros + 5 * n_unknowns
!       ia              2.2 * n_unknowns
!       u               2.2 * n_unknowns
!       f               2.2 * n_unknowns
!       kwork           5.4 * n_unknowns
!
!       BN: OK, that's how it once was, but I split the entire vector kwork
!       into three vectors: iw, ifg and icg, and pass these to subroutines
!       as work space.  (I already wrote this above.)
!
!    In the actual code (as recieved from K. Stueben), the values are bigger:
!
!     Vector         Needed length (guess)
!       a               6 * n_nonzeros + 12 * n_unknowns
!       ja              6 * n_nonzeros + 12 * n_unknowns
!       ia              3 * n_unknowns
!       u               3 * n_unknowns
!       f               3 * n_unknowns
!       kwork           6 * n_unknowns
!
!------------------------------------------------------------------------------!
!
!   Standard choices of parameters (as far as meaningful):
!
!     iswtch            =     4
!     iout              =     2
!     levelx            =    25
!     ifirst            =    13
!     ncyc              = 10110
!     eps               =     1.d-12
!     madapt            =    27
!     def_relax_down    =  1131
!     def_coarse_solver =   110
!     def_relax_up      =  1131
!     ecg1              =     0.
!     ecg2              =     0.25
!     ewt2              =     0.35
!     nwt               =     2
!     ntr               =     0
!
!   If any one of these parameters is 0 on input, its corresponding standard
!   value is used by Amg1r5.
!
!------------------------------------------------------------------------------!
!
!   Portability restrictions:
!
!   1. Most input parameters are composed of several digits, their
!      significance having been described above.  Be sure not to enter more
!      digits than your computer can store on an integer variable.
!      BN: I hope that, in the long run, I could use a different system.
!
!   2. Apart from Fortran intrinsic functions and service routines, there is
!      only one external reference to a program not contained in the Amg1r5
!      system, i.e. the linear system solver ndrv of the Yale Sparse Matrix
!      Package.  If you havn't access to this package, enter a dummy routine
!      ndrv and avoid choosing coarse_solver=2 (subparameter of def_coarse_solver).
!      Then ndrv isn't called by Amg1r5.  In this case, however, indefinite
!      problems will not be solvable.  The yale sparse matrix package is
!      freely available for non-profit purposes.  Contact the department of
!      computer science, Yale unitversity.  BN: I have removed the Yale SMP
!      package because it would not compile with modern Fortran compilers.
!
!   3. In Amg1r5 there is the parameter lratio, denoting the ratio of space
!      occupied by a double precision real variable and that of an integer.
!      For the IBM-version lratio has been set to 2.  Change this value if
!      necessary. (The same has to be done with the Yale SMP-routine ndrv.)
!      BN: OK, this is nowadays, when we have dynamic memory allocation and
!      when we don't re-use the same memory space by integers and reals,
!      rather obsolete.  Although I have introduced a parameter to store it
!      called AMG_L_RATIO, it is a candidate for deletion, really.
!
!------------------------------------------------------------------------------!
!
!   Authors:
!
!     John Ruge, Fort Collins (USA),
!       Institute for Computational Studies at CSU;
!
!     Klaus Stueben, D-5205 St. Augustin (W.-Germany),
!       Gesellschaft fuer Mathematik und Datenverarbeitung (GMD).
!
!     Rolf Hempel, D-5205 St. Augustin (W.-Germany),
!       Gesellschaft fuer Mathematik und Datenverarbeitung (GMD).
!
!   Re-factoring:
!
!     Bojan Niceno, CH-5232, Villigen-PSI (Switzerland)
!       Paul Scherrer Institute (PSI)
!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  real                    :: a_val(:)
  integer                 :: a_row(:), a_col(:)
  real                    :: phi(:), b(:)
  integer                 :: n
  real                    :: eps
  integer                 :: iswtch  ! Class 2 by Ruge-Stueben
!-----------------------------------[locals]-----------------------------------!
  integer :: digit(AMG_MAX_LEVELS)

  ! These arrays describe linear system
  real,    allocatable :: a(:), u(:), u_b(:), f(:), f_b(:)
  integer, allocatable :: ia(:)
  integer, allocatable :: ja(:)

  ! These were introduced later by B. Niceno to avoid the use of equivalence
  integer, allocatable :: iwork(:), jtr(:), iw(:), icg(:), ifg(:)

  integer :: n_nonzeros, n_unknowns, niw, nicg, nifg
  integer :: nda, ndu, ndw
  integer :: icgst, kevelx, kswtch, levels
  integer :: n_digits, ndicg
  integer :: levelx, ifirst, ncyc
  integer :: madapt, i,j,k,l, level
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !-------------------------------------------------------------!
  !   Process input argument iswtch (Class 2 by Ruge-Stueben)   !
  !-------------------------------------------------------------!
  if(iswtch .ne. 0) then
    kswtch = iswtch
  else
    kswtch = AMG_RUN_ALL_FOUR_STAGES
  end if

  !---------------------!
  !   Process eps too   !
  !---------------------!
  Amg % eps = eps

  !------------------------!
  !                        !
  !   If setup is needed   !
  !                        !
  !------------------------!
  if(kswtch .eq. AMG_RUN_ALL_FOUR_STAGES) then

    !   Work out the number of nonzeros and unknowns
    n_nonzeros = a_row(n+1)-1
    n_unknowns = n

    !----------------------------------------------------------!
    !   Initialize error to zero, assume there are no errors   !
    !----------------------------------------------------------!
    Amg % ierr = AMG_SUCCESS

    !-----------------------------------------------------------!
    !   Aditional (extended) dimensions needed for Amg solver   !
    !-----------------------------------------------------------!
    nda = 6 * n_nonzeros + 12 * n_unknowns
    ndu = 3 * n_unknowns
    ndw = 6 * n_unknowns

    ! This defines how (former) kwork is split
    icgst = n_unknowns+3
    ndicg = (ndw-icgst+1)/2
    if(ndicg .le. 0) then
      write(6, '(a)')  &
        ' *** error in Amg1r5: ndw too small ***'
      Amg % ierr = AMG_ERR_DIM_ICG_TOO_SMALL
      return
    end if
    niw  = icgst - 1
    nicg = ndicg
    nifg = ndw - icgst - ndicg + 1

    !----------------------------------------!
    !   Allocate memory for the AMG solver   !
    !----------------------------------------!
    if(.not. allocated(iwork)) then
      allocate(a(nda));      a(:)     = 0.0   ! matrix values
      allocate(ia(ndu));     ia(:)    = 0     ! row pointers
      allocate(ja(nda));     ja(:)    = 0     ! columns
      allocate(u(ndu));      u(:)     = 0.0
      allocate(u_b(ndu));    u_b(:)   = 0.0
      allocate(f(ndu));      f(:)     = 0.0
      allocate(f_b(ndu));    f_b(:)   = 0.0
      allocate(iw(niw));     iw(:)    = 0     ! used to be "ig"
      allocate(icg(nicg));   icg(:)   = 0
      allocate(ifg(nifg));   ifg(:)   = 0
      allocate(iwork(ndu));  iwork(:) = 0
      allocate(jtr(nda));    jtr(:)   = 0     ! has the same form as ja
    end if

    ! Copy the discretization to "Amg" workspace
    a (1:n_nonzeros)   = a_val(1:n_nonzeros)
    ia(1:n_unknowns+1) = a_row(1:n_unknowns+1)
    ja(1:n_nonzeros)   = a_col(1:n_nonzeros)
    u (1:n_unknowns)   = phi  (1:n_unknowns)
    f (1:n_unknowns)   = b    (1:n_unknowns)

    !--------------------------------------------------!
    !   Default values (Amg1r5: setup phase, output)   !
    !--------------------------------------------------!
    Amg % matrix = 12      ! rowsum /= 0.0; nonsymetric

    !--------------!
    !   Switches   !
    !--------------!
    ifirst       = 13        ! value from stuben: 13
    Amg % iout   =  1        ! can be from 0 to 4

    !-----------------------------------------------------------!
    !   More switches (these used to be in aux1r5 subroutine)   !
    !-----------------------------------------------------------!
    levelx                  = 0      ! just set kevelx to AMG_MAX_LEVELS
    ncyc                    = 12250  ! V, no CG, ||res|| < eps, 50 cycles
    madapt                  = 0
    Amg % def_relax_down    = 0      ! leave definition for later in Solve
    Amg % def_coarse_solver = AMG_SOLVER_CG
    Amg % def_relax_up      = 0      ! leave definition for later in Solve
    Amg % ecg1              = 0.0
    Amg % ecg2              = 0.25
    Amg % ewt2              = 0.35
    Amg % nwt               = 2
    Amg % ntr               = 0

    Amg % fine_solver = AMG_SOLVER_GS

    !------------------------------------!
    !   Limit kevelx to AMG_MAX_LEVELS   !
    !------------------------------------!
    if(levelx .gt. 0) then
      kevelx = min(levelx, AMG_MAX_LEVELS)
    else if (levelx .lt. 0) then
      write(6, '(a)')  &
        ' *** error in Amg1r5: illegal parameter levelx ***'
      Amg % ierr = AMG_ERR_LEVELX_INVALID
      return
    else
      kevelx = AMG_MAX_LEVELS
    end if

  !------------------------------------------!
  !                                          !
  !   Setup is not needed, but you have to   !
  !   copy the last solution and the right   !
  !   hand side to the local memory space    !
  !                                          !
  !------------------------------------------!
  else
    u(1:n_unknowns) = phi(1:n_unknowns)
    f(1:n_unknowns) = b  (1:n_unknowns)
  end if  ! kswtch

  !----------------------------------------------------------!
  !                                                          !
  !   Call whatever needs to be called according to kswtch   !
  !                                                          !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - +-------------!
  !   if kswtch is:                          call these stages:            !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !   AMG_RUN_ALL_FOUR_STAGES         (4):   Setup, First, Solve, Wrkcnt   !
  !   AMG_INITIALIZE_SOLVE_AND_REPORT (3):   -----, First, Solve, Wrkcnt   !
  !   AMG_SOLVE_AND_REPORT            (2):   -----, -----, Solve, Wrkcnt   !
  !   AMG_JUST_REPORT                 (1):   -----, -----, -----, Wrkcnt   !
  !------------------------------------------------------------------------!
  if(kswtch .eq. 4) then
    call Amg % Setup(n_unknowns, kevelx,  &
                     a, u, ia, ja,        &  ! linear system
                     iw, icg, ifg,        &  ! work arrays
                     levels,              &
                     iwork, jtr)

    Amg % imaxw(levels) = Amg % iminw(levels)
    !              1      81      81      81      81      81      81      8
    write(*, *)   ' level   min i   max i   size    min w   max w   total'
    do level = 1, levels
      write(*,'(99i8)')  &
        level, Amg % imin(level),  Amg % imax(level),       &
               Amg % imax(level) - Amg % imin(level) + 1,   &
               Amg % iminw(level), Amg % imaxw(level),      &
               Amg % imaxw(level) - Amg % iminw(level) + 1
    end do
    do level = 1, levels - 1
      do i = Amg % iminw(level), Amg % imaxw(level)
        k = iw(i)  ! position in the matrix
        j = ja(k)  ! column number
        l = Amg % Level_Of_Cell(j)
        ! write(*, '(99i8)')  k, j, l
        Assert(l == level+1)
      end do
    end do

    !------------------------------------!
    !   Reduce excessive size matrices   !
    !- - - - - - - - - - - - - - - - - - +----------------------------!
    !   Although this seems cool, it presumes that only one call to   !
    !   Setup will be made, with subsequent calls skipping it. I am   !
    !   not sure how to work around it. maybe a more complex logic    !
    !   will be needed at the section for memory allocation above?    !
    !-----------------------------------------------------------------!

    do level = 1, levels
      call Amg % Store_Level(level, levels, a, u, f, ia, ja, iw, icg, ifg)
    end do

    write(6, '(a)')  'De-allocating iwork and jtr, iw, icg and ifg'
    deallocate(iwork, jtr, iw, icg, ifg)

    write(6, '(a)')  'De-allocating a, ia and ja'
    deallocate(a, ia, ja)

    write(6, '(a,i11,a,i11)')  &
      'Reducing u   from ', size(u, 1), ' to ', Amg % imax(levels)
    call Amg % Reduce_Real(u, Amg % imax(levels))

    write(6, '(a,i11,a,i11)')  &
      'Reducing f   from ', size(f, 1), ' to ', Amg % imax(levels)
    call Amg % Reduce_Real(f, Amg % imax(levels))

    if(Amg % ierr .gt. 0) return
  end if

  if(kswtch .ge. 3) then
    call Amg % First_Guess(ifirst, u)
  end if

  if(kswtch .ge. 2) then
    call Amg % Update_U_And_F_At_Level(1, vec_u=u, vec_f=f, for_real=.true.)
    call Amg % Solve(madapt, ncyc,               &
                     a, u, u_b, f, f_b, ia, ja,  &  ! linear system
                     iw, icg, ifg,               &  ! work arrays
                     levels)
    if(Amg % ierr .gt. 0) return
  end if

  if(kswtch .ge. 1) then
    call Amg % Wrkcnt(ia, iw, levels)
  end if

  !------------------------!
  !   Fetch the solution   !
  !------------------------!
  phi(1:n_unknowns) = Amg % lev(1) % u(1:n_unknowns)

  if(kswtch .lt. 1 .or. kswtch .gt. 4) then
    write(6, '(a)')  &
      ' *** error in Amg1r5: illegal parameter Amg % iswtch ***'
    Amg % ierr = AMG_ERR_ISWTCH_INVALID
  end if

  end subroutine
