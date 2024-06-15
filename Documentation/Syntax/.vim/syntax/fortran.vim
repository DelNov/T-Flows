"=====================================================================[T-Flows]=
" The idea here is as follows:
"
" Expand the existing fortran syntax file which comes with vim, with code
" contructs pertinent to T-Flow.  Ideally, that should be followed by
" setting or selecting a colorsheme in vim, but defining a color scheme
" seems to be a daunting task.
"
" Rather, fine-tuning of colors in vim is don through .vimrc file, using
" the command highlight (typical format :highlight String ctermfg=Green)
"
" There is a CATCH with setting highlights in .vimrc, though: formatting
" elements in .vimrc should NOT be set in this file!!!  (Not sure anymore)
"
" A sample snippet from .vimrc which worked well now follows:
"
" from .vimrc: "----------------------------------------------------------------
" from .vimrc: " With the following commands, I am trying to tune color scheme
" from .vimrc: "
" from .vimrc: " Commands which might prove to be useful in the future
" from .vimrc: "
" from .vimrc: " Preprocessor directives
" from .vimrc: " :highlight PreProc  ctermfg=Red
" from .vimrc: "
" from .vimrc: " Everything else:
" from .vimrc: " :highlight Normal  ctermfg=Yellow
" from .vimrc: "----------------------------------------------------------------
" from .vimrc: 
" from .vimrc: " Keyords should be Fortran commands
" from .vimrc: highlight fortranKeyword ctermfg=Yellow
" from .vimrc: 
" from .vimrc: " Ghost numbers are not desirable, mark them red
" from .vimrc: highlight fortranNumber ctermfg=Red
" from .vimrc: highlight fortranFloat  ctermfg=Red
" from .vimrc: 
" from .vimrc: " Constants are better than ghost numbers, mark them blue ...
" from .vimrc: " ... they are safer, they are in effect former ghost numbers
" from .vimrc: highlight fortranConstant ctermfg=Blue
" from .vimrc: 
" from .vimrc: " Continuation lines.  I like them visible and a bit psychedelic
" from .vimrc: highlight fortranContinueMark ctermfg=Magenta
" from .vimrc: 
" from .vimrc: highlight fortranString ctermfg=Gray
" from .vimrc: 
" from .vimrc: " Operators ... I don't know, maybe the same as keywords
" from .vimrc: highlight fortranOperator ctermfg=Yellow
" from .vimrc: 
" from .vimrc: " Just Comment didn't work, fortranComment did
" from .vimrc: highlight fortranComment ctermfg=LightBlue
" from .vimrc: 
" from .vimrc: " Try to over-ride Fortran objects in T-Flows
" from .vimrc: highlight fortranObjectTflows cterm=Bold ctermfg=Green
" from .vimrc: 
" from .vimrc: " MPI calls in T-Flows - they are few and far in between, yet dangerous
" from .vimrc: highlight fortranMpiTflows ctermfg=Red
"---------------------------------------------------------------------[T-Flows]-

" Vim syntax file
" Language:	Fortran 2008 (and older: Fortran 2003, 95, 90, and 77)
" Version:	(v104) 2021 April 06
" Maintainer:	Ajit J. Thakkar <ajit@unb.ca>; <http://www2.unb.ca/~ajit/>
" Usage:	For instructions, do :help fortran-syntax from Vim
" Credits:
"  Version 0.1 for Fortran 95 was created in April 2000 by Ajit Thakkar from an
"  older Fortran 77 syntax file by Mario Eusebio and Preben Guldberg.
"  Since then, useful suggestions and contributions have been made, in order, by:
"  Andrej Panjkov, Bram Moolenaar, Thomas Olsen, Michael Sternberg, Christian Reile,
"  Walter Dieudonne, Alexander Wagner, Roman Bertle, Charles Rendleman,
"  Andrew Griffiths, Joe Krahn, Hendrik Merx, Matt Thompson, Jan Hermann,
"  Stefano Zaghi, Vishnu V. Krishnan, Judicael Grasset, Takuma Yoshida,
"  Eisuke Kawashima, Andre Chalella, and Fritz Reese.

if exists("b:current_syntax")
  finish
endif

let s:cpo_save = &cpo
set cpo&vim

" Choose fortran_dialect using the priority:
" source file directive > buffer-local value > global value > file extension
" first try using directive in first three lines of file
let b:fortran_retype = getline(1)." ".getline(2)." ".getline(3)
if b:fortran_retype =~? '\<fortran_dialect\s*=\s*F\>'
  let b:fortran_dialect = "F"
elseif b:fortran_retype =~? '\<fortran_dialect\s*=\s*f08\>'
  let b:fortran_dialect = "f08"
elseif !exists("b:fortran_dialect")
  if exists("g:fortran_dialect") && g:fortran_dialect =~# '\<F\|f08\>'
    " try global variable
    let b:fortran_dialect = g:fortran_dialect
  else         " nothing found, so use default
    let b:fortran_dialect = "f08"
  endif
endif
unlet! b:fortran_retype
" make sure buffer-local value is not invalid
if b:fortran_dialect !~# '\<F\|f08\>'
  let b:fortran_dialect = "f08"
endif

"=====================================================================[T-Flows]=
" Do not ignore the case
" syn case ignore
"---------------------------------------------------------------------[T-Flows]-

syn match fortranConstructName   "^\s*\zs\a\w*\ze\s*:"
if exists("fortran_more_precise")
  syn match fortranConstructName   "\(\<end\s*do\s\+\)\@11<=\a\w*"
  syn match fortranConstructName   "\(\<end\s*if\s\+\)\@11<=\a\w*"
  syn match fortranConstructName   "\(\<end\s*select\s\+\)\@15<=\a\w*"
endif

syn match fortranUnitHeader        "\<end\>"
syn match fortranType              "\<character\>"
syn match fortranType              "\<complex\>"
syn match fortranType              "\<integer\>"
syn match fortranType              "\<real\>"
syn match fortranType              "\<logical\>"
syn keyword fortranType            intrinsic
syn match fortranType              "\<implicit\>"
syn keyword fortranStructure       dimension
syn keyword fortranStorageClass    parameter save
syn match fortranUnitHeader        "\<subroutine\>"
"=====================================================================[T-Flows]=
" Got rid of the special group for call ...
" syn keyword fortranCall		call
"---------------------------------------------------------------------[T-Flows]-
syn match fortranUnitHeader  "\<function\>"
syn match fortranUnitHeader  "\<program\>"
syn match fortranUnitHeader  "\<block\>"
"=====================================================================[T-Flows]=
" ... and added it (call) to the keywords
" syn keyword fortranKeyword	return stop
syn keyword fortranKeyword	return stop call
"---------------------------------------------------------------------[T-Flows]-
syn keyword fortranConditional     else then
syn match fortranConditional       "\<if\>"
syn match fortranConditionalOb     "\<if\s*(.*)\s*\d\+\s*,\s*\d\+\s*,\s*\d\+\s*$"
syn match fortranRepeat            "\<do\>"

syn keyword fortranTodo		contained todo fixme

"Catch errors caused by too many right parentheses
syn region fortranParen transparent start="(" end=")" contains=ALLBUT,fortranParenError,@fortranCommentGroup,cIncluded,@spell
syn match  fortranParenError   ")"

syn match fortranOperator	"\.\s*n\=eqv\s*\."
syn match fortranOperator	"\.\s*\(and\|or\|not\)\s*\."
syn match fortranOperator	"\(+\|-\|/\|\*\)"
syn match fortranTypeOb		"\<character\s*\*"

syn match fortranBoolean	"\.\s*\(true\|false\)\s*\."

syn keyword fortranReadWrite	backspace close endfile inquire open print read rewind write

"If tabs are allowed then the left margin checks do not work
if exists("fortran_have_tabs")
  syn match fortranTab		"\t"  transparent
else
  syn match fortranTab		"\t"
endif

"=====================================================================[T-Flows]=
" Removed name and file from the list below, too often used as variables
"---------------------------------------------------------------------[T-Flows]-
syn keyword fortranIO  access blank direct exist fmt form formatted iostat named
syn keyword fortranIO  nextrec number opened rec recl sequential status unformatted unit

syn keyword fortranIntrinsicR		alog alog10 amax0 amax1 amin0 amin1 amod cabs ccos cexp clog csin csqrt dabs dacos dasin datan datan2 dcos dcosh ddim dexp dint dlog dlog10 dmax1 dmin1 dmod dnint dsign dsin dsinh dsqrt dtan dtanh float iabs idim idint idnint ifix isign max0 max1 min0 min1 sngl

" Intrinsics provided by some vendors
syn keyword fortranExtraIntrinsic	algama cdabs cdcos cdexp cdlog cdsin cdsqrt cqabs cqcos cqexp cqlog cqsin cqsqrt dcmplx dconjg derf derfc dfloat dgamma dimag dlgama iqint qabs qacos qasin qatan qatan2 qcmplx qconjg qcos qcosh qdim qerf qerfc qexp qgamma qimag qlgama qlog qlog10 qmax1 qmin1 qmod qnint qsign qsin qsinh qsqrt qtan qtanh

syn keyword fortranKeyword   abs acos aimag aint anint asin atan atan2 char cmplx conjg cos cosh exp
syn keyword fortranKeyword   ichar index int log log10 max min nint sign sin sinh sqrt tan tanh
"=====================================================================[T-Flows]=
" Changed cast functions from Intrinsic to Fortran keywords
"---------------------------------------------------------------------[T-Flows]-
syn match fortranKeyword     "\<len\s*[(,]"me=s+3
syn match fortranKeyword     "\<real\s*("me=s+4
syn match fortranKeyword     "\<logical\s*("me=s+7
syn match fortranType        "\<implicit\s\+real\>"
syn match fortranType        "\<implicit\s\+logical\>"

"Numbers of various sorts
" Integers
syn match fortranNumber	display "\<\d\+\(_\a\w*\)\=\>"
" floating point number, without a decimal point
syn match fortranFloatIll	display	"\<\d\+[deq][-+]\=\d\+\(_\a\w*\)\=\>"
" floating point number, starting with a decimal point
syn match fortranFloatIll	display	"\.\d\+\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, no digits after decimal
syn match fortranFloatIll	display	"\<\d\+\.\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, D or Q exponents
syn match fortranFloatIll	display	"\<\d\+\.\d\+\([dq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number
syn match fortranFloat	display	"\<\d\+\.\d\+\(e[-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" binary number
syn match fortranBinary	display	"b["'][01]\+["']"
" octal number
syn match fortranOctal	display	"o["'][0-7]\+["']"
" hexadecimal number
syn match fortranHex	display	"z["'][0-9A-F]\+["']"
" Numbers in formats
syn match fortranFormatSpec	display	"\d*f\d\+\.\d\+"
syn match fortranFormatSpec	display	"\d*e[sn]\=\d\+\.\d\+\(e\d+\>\)\="
syn match fortranFormatSpec	display	"\d*\(d\|q\|g\)\d\+\.\d\+\(e\d+\)\="
syn match fortranFormatSpec	display	"\d\+x\>"
" The next match cannot be used because it would pick up identifiers as well
" syn match fortranFormatSpec	display	"\<\(a\|i\)\d\+"

" Numbers as labels
syn match fortranLabelNumber	display	"^\d\{1,5}\s"me=e-1
syn match fortranLabelNumber	display	"^ \d\{1,4}\s"ms=s+1,me=e-1
syn match fortranLabelNumber	display	"^  \d\{1,3}\s"ms=s+2,me=e-1
syn match fortranLabelNumber	display	"^   \d\d\=\s"ms=s+3,me=e-1
syn match fortranLabelNumber	display	"^    \d\s"ms=s+4,me=e-1

if exists("fortran_more_precise")
  " Numbers as targets
  syn match fortranTarget	display	"\(\<if\s*(.\+)\s*\)\@<=\(\d\+\s*,\s*\)\{2}\d\+\>"
  syn match fortranTarget	display	"\(\<do\s\+\)\@11<=\d\+\>"
  syn match fortranTarget	display	"\(\<go\s*to\s*(\=\)\@11<=\(\d\+\s*,\s*\)*\d\+\>"
endif

syn keyword fortranTypeR           external
syn keyword fortranIOR             format
syn match fortranKeywordR          "\<continue\>"
syn match fortranKeyword           "^\s*\d\+\s\+continue\>"
syn match fortranKeyword           "\<go\s*to\>"
syn match fortranKeywordDel        "\<go\s*to\ze\s\+.*,\s*(.*$"
syn match fortranKeywordOb         "\<go\s*to\ze\s*(\d\+.*$"
syn region fortranStringR          start=+'+ end=+'+ contains=fortranContinueMark,fortranLeftMargin,fortranSerialNumber
syn keyword fortranIntrinsicR      dim lge lgt lle llt mod
syn keyword fortranKeywordDel      assign pause

syn match fortranType        "\<type\>"
syn keyword fortranType      none

syn keyword fortranStructure       private public intent optional
syn keyword fortranStructure       pointer target allocatable
syn keyword fortranStorageClass    in out
syn match fortranStorageClass      "\<kind\s*="me=s+4
syn match fortranStorageClass      "\<len\s*="me=s+3

syn match fortranUnitHeader        "\<module\>"
syn match fortranUnitHeader        "\<submodule\>"
syn keyword fortranUnitHeader      use only contains
syn keyword fortranUnitHeader      result operator assignment
syn match fortranUnitHeader        "\<interface\>"
syn keyword fortranKeyword         allocate deallocate nullify cycle exit
syn match fortranConditional       "\<select\>"
syn keyword fortranConditional     case default where elsewhere

syn match fortranOperator    "\(\(>\|<\)=\=\|==\|/=\|=\)"
syn match fortranOperator    "=>"

syn region fortranString     start=+"+ end=+"+	contains=fortranLeftMargin,fortranContinueMark,fortranSerialNumber
syn keyword fortranIO        pad position action delim readwrite
syn keyword fortranIO        eor advance nml

syn keyword fortranKeyword   adjustl adjustr all allocated any associated bit_size btest
syn keyword fortranKeyword   ceiling count cshift date_and_time digits dot_product eoshift
syn keyword fortranKeyword   epsilon exponent floor fraction huge iand ibclr ibits ibset
syn keyword fortranKeyword   ieor ior ishft ishftc lbound len_trim matmul maxexponent
syn keyword fortranKeyword   maxloc maxval merge minexponent minloc minval modulo mvbits
syn keyword fortranKeyword   nearest pack precision present product radix random_number
syn keyword fortranKeyword   random_seed range repeat reshape rrspacing
"=====================================================================[T-Flows]=
" Changed a few more functions from Intrinsic to Fortran keywords
"---------------------------------------------------------------------[T-Flows]-
syn keyword fortranKeyword   scale scan selected_int_kind selected_real_kind set_exponent
syn keyword fortranKeyword   shape size spacing spread sum system_clock tiny transpose
syn keyword fortranKeyword   trim ubound unpack verify
syn match fortranIntrinsic   "\<not\>\(\s*\.\)\@!"me=s+3
syn match fortranIntrinsic   "\<kind\>\s*[(,]"me=s+4

syn match  fortranUnitHeader       "\<end\s*function"
syn match  fortranUnitHeader       "\<end\s*interface"
syn match  fortranUnitHeader       "\<end\s*module"
syn match  fortranUnitHeader       "\<end\s*submodule"
syn match  fortranUnitHeader       "\<end\s*program"
syn match  fortranUnitHeader       "\<end\s*subroutine"
syn match  fortranUnitHeader       "\<end\s*block"
syn match  fortranRepeat           "\<end\s*do"
syn match  fortranConditional      "\<end\s*where"
syn match  fortranConditional      "\<select\s*case"
syn match  fortranConditional      "\<end\s*select"
syn match  fortranType             "\<end\s*type"
syn match  fortranType             "\<in\s*out"

syn keyword fortranType            procedure
syn match  fortranType             "\<module\ze\s\+procedure\>"
syn keyword fortranIOR             namelist
syn keyword fortranConditionalR    while
syn keyword fortranIntrinsicR      achar iachar transfer

syn keyword fortranInclude         include
syn keyword fortranStorageClassR   sequence

syn match   fortranConditional     "\<end\s*if"
syn match   fortranIO              contains=fortranOperator "\<e\(nd\|rr\)\s*=\s*\d\+"
syn match   fortranConditional     "\<else\s*if"

syn keyword fortranUnitHeaderOb    entry
syn match fortranTypeR             display "double\s\+precision"
syn match fortranTypeR             display "double\s\+complex"
syn match fortranUnitHeaderR       display "block\s\+data"
syn keyword fortranStorageClassR   common equivalence data
syn keyword fortranIntrinsicR      dble dprod
syn match   fortranOperatorR       "\.\s*[gl][et]\s*\."
syn match   fortranOperatorR       "\.\s*\(eq\|ne\)\s*\."

syn keyword fortranRepeat          forall
syn match fortranRepeat            "\<end\s*forall"
syn keyword fortranIntrinsic       null cpu_time
syn match fortranType              "\<elemental\>"
syn match fortranType              "\<pure\>"
syn match fortranType              "\<impure\>"
syn match fortranType              "\<recursive\>"
if exists("fortran_more_precise")
  syn match fortranConstructName   "\(\<end\s*forall\s\+\)\@15<=\a\w*\>"
endif

if b:fortran_dialect == "f08"
  " F2003
  syn keyword fortranKeyword       command_argument_count get_command get_command_argument get_environment_variable
  syn keyword fortranKeyword       is_iostat_end is_iostat_eor move_alloc new_line selected_char_kind same_type_as
  syn keyword fortranKeyword       extends_type_of
  " ISO_C_binding
  syn keyword fortranConstant      c_null_char c_alert c_backspace c_form_feed c_new_line c_carriage_return
  syn keyword fortranConstant      c_horizontal_tab c_vertical_tab c_int c_short c_long c_long_long c_signed_char
  syn keyword fortranConstant      c_size_t c_int8_t c_int16_t c_int32_t c_int64_t c_int_least8_t c_int_least16_t
  syn keyword fortranConstant      c_int_least32_t c_int_least64_t c_int_fast8_t c_int_fast16_t c_int_fast32_t
  syn keyword fortranConstant      c_int_fast64_t c_intmax_t C_intptr_t c_float c_double c_long_double c_float_complex
  syn keyword fortranConstatn      c_double_complex c_long_double_complex c_bool c_char c_null_ptr c_null_funptr
  syn keyword fortranIntrinsic     iso_c_binding c_loc c_funloc c_associated  c_f_pointer c_f_procpointer
  syn keyword fortranType          c_ptr c_funptr
  " ISO_Fortran_env
  syn keyword fortranConstant      iso_fortran_env character_storage_size error_unit file_storage_size input_unit iostat_end iostat_eor numeric_storage_size output_unit
  " IEEE_arithmetic
  syn keyword fortranIntrinsic     ieee_arithmetic ieee_support_underflow_control ieee_get_underflow_mode ieee_set_underflow_mode

  syn keyword fortranReadWrite     flush wait
  syn keyword fortranIO            decimal round iomsg
  syn keyword fortranType          asynchronous nopass non_overridable pass protected volatile extends import
  syn keyword fortranType          non_intrinsic value bind deferred generic final enumerator
  syn match fortranType            "\<abstract\>"
  syn match fortranType            "\<class\>"
  syn match fortranType            "\<associate\>"
  syn match fortranType            "\<end\s*associate"
  syn match fortranType            "\<enum\s*,\s*bind\s*(\s*c\s*)"
  syn match fortranType            "\<end\s*enum"
  syn match fortranConditional     "\<select\s*type"
  syn match fortranConditional     "\<type\s*is\>"
  syn match fortranConditional     "\<class\s*is\>"
  syn match fortranUnitHeader      "\<abstract\s*interface\>"
  syn match fortranOperator        "\([\|]\)"

  " F2008
  syn keyword fortranKeyword       acosh asinh atanh bessel_j0 bessel_j1 bessel_jn bessel_y0 bessel_y1 bessel_yn
  syn keyword fortranKeyword       erf erfc erfc_scaled gamma log_gamma hypot norm2
  syn keyword fortranKeyword       atomic_define atomic_ref execute_command_line leadz trailz storage_size merge_bits
  syn keyword fortranKeyword       bge bgt ble blt dshiftl dshiftr findloc iall iany iparity image_index lcobound ucobound maskl maskr num_images parity popcnt poppar shifta shiftl shiftr this_image
  syn keyword fortranIO            newunit
  syn keyword fortranType          contiguous
  syn keyword fortranRepeat        concurrent

"=====================================================================[T-Flows]=
" T-Flows specific
" Start with constants from Const_Mod:
  syn keyword fortranConstant      PROGRAM_NAME
  syn keyword fortranConstant      VL  SL  DL   DP  SP  IP  LP  RP  MAX_STRING_ITEMS  MAX_PETSC_MEMBERS
  syn keyword fortranConstant      VERSION_CFN  VERSION_DIM  VERSION_BACKUP
  syn keyword fortranConstant      YOCTO  ZEPTO  ATTO  FEMTO  PICO  NANO  MICRO  MILI
  syn keyword fortranConstant      YOTTA  ZETTA  EXA   PETA   TERA  GIGA  MEGA   KILO
  syn keyword fortranConstant      HUGE  TINY  HUGE_INT  EULER  PI
  syn keyword fortranConstant      ONE_THIRD  TWO_THIRDS  ONE_SIXTH
  syn keyword fortranConstant      MD  MAX_VARS_INTERFACE  MAX_TURBULENT_PLANES
  syn keyword fortranConstant      GROWTH_MARGIN
" Constants from Tokenizer_Mod
  syn keyword fortranConstant      MAX_TOKENS
" Constants related to indentation (I think they are defined in a few places - bad!
  syn keyword fortranConstant      IN_0  IN_1  IN_2  IN_3  IN_4  IN_5
" Constants from Info_Mod
  syn keyword fortranConstant      L_LINE  L_BOX  MAX_USER_LINES
" Constants from Numerics_Mod
  syn keyword fortranConstant      UPWIND  CENTRAL  LUDS  QUICK  SMART  GAMMA
  syn keyword fortranConstant      MINMOD  BLENDED  SUPERBEE  AVL_SMART  CICSAM  STACS
  syn keyword fortranConstant      LINEAR  PARABOLIC  RUNGE_KUTTA_3
  syn keyword fortranConstant      SIMPLE  PISO  CHOI
  syn keyword fortranConstant      LEAST_SQUARES  GAUSS_THEOREM
" Constants from Turb_Mod
  syn keyword fortranConstant      NO_TURBULENCE_MODEL  DNS  LES_SMAGORINSKY  LES_DYNAMIC
  syn keyword fortranConstant      LES_WALE  LES_TVM  K_EPS  K_EPS_ZETA_F  DES_SPALART
  syn keyword fortranConstant      SPALART_ALLMARAS  RSM_HANJALIC_JAKIRLIC  RSM_MANCEAU_HANJALIC
  syn keyword fortranConstant      HYBRID_LES_RANS HYBRID_LES_PRANDTL  STABILIZED  SGDH  GGDH  AFM  A_POW  B_POW
  syn keyword fortranConstant      SWITCH_DISTANCE  SWITCH_VELOCITY  THERMALLY_DRIVEN  DENSITY_DRIVEN  NO_BUOYANCY
" Constants from Swarm_Mod
  syn keyword fortranConstant      BROWNIAN_FUKAGATA  DISCRETE_RANDOM_WALK  N_I_VARS  N_L_VARS  N_R_VARS
" Constants from Region_Mod
  syn keyword fortranConstant      INFLOW  WALL  WALLFL  OUTFLOW  SYMMETRY  CONVECT  PRESSURE
  syn keyword fortranConstant      INSIDE  BUFFER  PERIODIC_X  PERIODIC_Y  PERIODIC_Z  UNDEFINED
" Constants from Solver_Mod and other related modules
  syn keyword fortranConstant      NATIVE  PETSC  PETSC_ACTIVE  OUT_OF_ITS
" Constants from various procedures and modules
  syn keyword fortranConstant      YES  NO  DEBUG  BEGIN  DEFAULT_TOLERANCE
" Constants from File_Mod
  syn keyword fortranConstant      BUFFER_SIZE  MAX_ITEMS
" Constant from Grid_Mod
  syn keyword fortranConstant      MAX_CLUSTERS
" Constant from Convert_Mod
  syn keyword fortranConstant      MAX_N  ITERS  FOR_BUILDINGS  FOR_POROSITIES  FOR_CHIMNEYS
" Constant from Native_Mod
  syn keyword fortranConstant      MATRIX_ONE  MATRIX_UVW  MATRIX_PP   MATRIX_T
" After the constants, I have alternating definitions of types and objects derived from them
  syn keyword fortranTypeTflows    Domain_Type    Point_Type     Block_Type     Line_Type         Range_Type     Read_Controls_Type
  syn keyword fortranObjectTflows  Dom            points Point   blocks         lines             ranges         Read_Control  Rc
  syn keyword fortranTypeTflows    Generate_Type  Convert_Type   Divide_Type    Grid_Type         Control_Type   Time_Type
  syn keyword fortranObjectTflows  Generate       Convert        Divide         Grid  Prim  Dual  Control        Time
  syn keyword fortranTypeTflows    Math_Type      Sort_Type      File_Type      String_Type       Work_Type      Tokenizer_Type
  syn keyword fortranObjectTflows  Math           Sort           File           String            Work           Line  Tok
  syn keyword fortranTypeTflows    Comm_Type      Backup_Type    Field_Type     Turb_Type         Vof_Type       Swarm_Type
  syn keyword fortranObjectTflows  Comm  Global   Backup  Bac    Flow  Fld      Turb  Tur         Vof            Swarm Swr
  syn keyword fortranTypeTflows    Bulk_Type      Face_Type      Iter_Type      Metis_Type        Omp_Type
  syn keyword fortranObjectTflows                 v_flux         Iter           Metis             Omp
  syn keyword fortranTypeTflows    Front_Type     Surf_Type      Elem_Type      Side_Type         Vert_Type      Particle_Type
  syn keyword fortranObjectTflows  Front          Surf           Elem           side              Vert           Part  Particle
  syn keyword fortranTypeTflows    Monitor_Type   Results_Type   Porosity_Type  Profiler_Type     Message_Type   Info_Type
  syn keyword fortranObjectTflows  Monitor        Results        Por            Profiler  Prof    Message  Msg   Info
  syn keyword fortranTypeTflows    Matrix_Type    Vector_Type    Solver_Type    Native_Type       Petsc_Type     Work_Petsc_Type
  syn keyword fortranObjectTflows  A M C Mat      vector         Sol            Nat               Pet            Work_Pet  Member
  syn keyword fortranTypeTflows    Process_Type   Pattern_Type   Isoap_Type     Polyhedron_Type   Stl_Type       Iso_Polygons_Type
  syn keyword fortranObjectTflows  Process        Pat            Isoap          Polyhedron  Pol   Stl            Iso_Polygons  Iso
  syn keyword fortranTypeTflows    Gpu_Type
  syn keyword fortranObjectTflows  Gpu
  syn keyword fortranTypeTflows    Var_Type
  syn keyword fortranObjectTflows  u  v  w  ui  uj  uk  p  pp  t  kin  eps  zeta  f22  uu  vv  ww  uv  vw  uw  ut  vt  wt  t2  vis  phi
  syn keyword fortranTypeTflows    Eddy_Type  Eddies_Type  Turb_Plane_Type  Memory_Type
  syn keyword fortranObjectTflows  Eddies  Plane  Turb_Planes               Enlarge  Mem
" Items which follow are not really objects, but I don't know where else to put them
  syn keyword fortranObjectTflows  This_Proc  N_Procs  First_Proc  Sequential_Run  Parallel_Run
  syn keyword fortranPetscTflows  C_Petsc_Finalize  C_Petsc_Initialize  C_Petsc_Ksp_Converged_Reason  C_Petsc_Ksp_Create  C_Petsc_Ksp_Destroy
  syn keyword fortranPetscTflows  C_Petsc_Ksp_Get_Iteration_Number  C_Petsc_Ksp_Get_Residual_Norm C_Petsc_Ksp_Set_Initial_Guess_Nonzero
  syn keyword fortranPetscTflows  C_Petsc_Ksp_Set_Operators  C_Petsc_Ksp_Set_Preconditioner  C_Petsc_Ksp_Set_Tolerances  C_Petsc_Ksp_Set_Type
  syn keyword fortranPetscTflows  C_Petsc_Ksp_Solve  C_Petsc_Log_Default_Begin  C_Petsc_Log_View  C_Petsc_Mat_Aij_Set_Preallocation
  syn keyword fortranPetscTflows  C_Petsc_Mat_Assemble  C_Petsc_Mat_Create  C_Petsc_Mat_Destroy  C_Petsc_Mat_Remove_Null_Space
  syn keyword fortranPetscTflows  C_Petsc_Mat_Set_Null_Space  C_Petsc_Mat_Set_Sizes  C_Petsc_Mat_Set_Type_To_Mat_Aij  C_Petsc_Mat_Set_Value
  syn keyword fortranPetscTflows  C_Petsc_Options_Set_Value  C_Petsc_Vec_Assemble  C_Petsc_Vec_Create_Mpi  C_Petsc_Vec_Destroy  C_Petsc_Vec_Get_Values
  syn keyword fortranPetscTflows  C_Petsc_Vec_Set_Value  C_Petsc_Mat_Zero_Entries
" Macros might need a special group of their own.
  syn keyword fortranMacroTflows   Boundary_Regions  Boundary_And_Inside_Regions  Boundary_Inside_And_Buffer_Regions
  syn keyword fortranMacroTflows   All_Regions  Faces_In_Region  Faces_In_Domain_And_At_Buffers
  syn keyword fortranMacroTflows   Cells_In_Region  Cells_In_Domain  Cells_In_Domain_And_Buffers  Cells_In_Buffers  Cell_In_This_Proc
  syn keyword fortranMacroTflows   Cells_At_Boundaries  Cells_At_Boundaries_In_Domain_And_Buffers
  syn keyword fortranMacroTflows   Assert  Unused
" Finally, a few global functions which I don't really like in the code
  syn keyword fortranGlobalTflows  Adjust_Dim  Adjust_First_Dim  Swap_Int  Swap_Real  Key_Ind
"---------------------------------------------------------------------[T-Flows]-

"==============================================================[MPI in T-Flows]=
" Here are MPI calls from T-Flows
" Note that they are not set in this file at all, only in the .vimrc
  syn keyword fortranMpiTflows     Mpi_Init  Mpi_Comm_Size  Mpi_Comm_Rank  Mpi_Barrier  Mpi_Finalize
  syn keyword fortranMpiTflows     Mpi_Allreduce  Mpi_File_Close  Mpi_File  Mpi_File_Set_View Mpi_Sendrecv_Replace
  syn keyword fortranMpiTflows     Mpi_Send  Mpi_Sendrecv  Mpi_Status  Mpi_Recv  Mpi_File_Write  Mpi_Datatype
  syn keyword fortranMpiTflows     Mpi_Type_Create_Indexed_Block  Mpi_Type_Commit
  syn keyword fortranMpiTflows     Mpi_File_Open  Mpi_File_Read
  syn keyword fortranMpiTflows     MPI_COMM_WORLD  MPI_CHARACTER  MPI_INFO_NULL  MPI_STATUS_IGNORE  MPI_IN_PLACE
  syn keyword fortranMpiTflows     MPI_INTEGER  MPI_LOGICAL  MPI_DOUBLE_PRECISION  MPI_REAL
  syn keyword fortranMpiTflows     MPI_LOR  MPI_SUM  MPI_MAX  MPI_MIN  MPI_MODE_WRONLY  MPI_MODE_CREATE
"--------------------------------------------------------------[MPI in T-Flows]-

" CUDA fortran
  syn keyword fortranTypeCUDA      host global device value
  syn keyword fortranTypeCUDA      shared constant pinned texture
  syn keyword fortranTypeCUDA      dim1 dim2 dim3 dim4
  syn keyword fortranTypeCUDA      cudadeviceprop cuda_count_kind cuda_stream_kind
  syn keyword fortranTypeCUDA      cudaEvent cudaFuncAttributes cudaArrayPtr
  syn keyword fortranTypeCUDA      cudaSymbol cudaChannelFormatDesc cudaPitchedPtr
  syn keyword fortranTypeCUDA      cudaExtent cudaMemcpy3DParms
  syn keyword fortranTypeCUDA      cudaFuncCachePreferNone cudaFuncCachePreferShared
  syn keyword fortranTypeCUDA      cudaFuncCachePreferL1 cudaLimitStackSize
  syn keyword fortranTypeCUDA      cudaLimitPrintfSize cudaLimitMallocHeapSize
  syn keyword fortranTypeCUDA      cudaSharedMemBankSizeDefault cudaSharedMemBankSizeFourByte cudaSharedMemBankSizeEightByte
  syn keyword fortranTypeCUDA      cudaEventDefault cudaEventBlockingSync cudaEventDisableTiming
  syn keyword fortranTypeCUDA      cudaMemcpyHostToDevice cudaMemcpyDeviceToHost
  syn keyword fortranTypeCUDA      cudaMemcpyDeviceToDevice
  syn keyword fortranTypeCUDA      cudaErrorNotReady cudaSuccess cudaErrorInvalidValue
  syn keyword fortranTypeCUDA      c_devptr

  syn match fortranStringCUDA      "blockidx%[xyz]"
  syn match fortranStringCUDA      "blockdim%[xyz]"
  syn match fortranStringCUDA      "griddim%[xyz]"
  syn match fortranStringCUDA      "threadidx%[xyz]"

  syn keyword fortranIntrinsicCUDA    warpsize syncthreads syncthreads_and syncthreads_count syncthreads_or threadfence threadfence_block threadfence_system gpu_time allthreads anythread ballot
  syn keyword fortranIntrinsicCUDA    atomicadd atomicsub atomicmax atomicmin atomicand atomicor atomicxor atomicexch atomicinc atomicdec atomiccas sizeof __shfl __shfl_up __shfl_down __shfl_xor
  syn keyword fortranIntrinsicCUDA    cudaChooseDevice cudaDeviceGetCacheConfig cudaDeviceGetLimit cudaDeviceGetSharedMemConfig cudaDeviceReset cudaDeviceSetCacheConfig cudaDeviceSetLimit cudaDeviceSetSharedMemConfig cudaDeviceSynchronize cudaGetDevice cudaGetDeviceCount cudaGetDeviceProperties cudaSetDevice cudaSetDeviceFlags cudaSetValidDevices
  syn keyword fortranIntrinsicCUDA    cudaThreadExit cudaThreadSynchronize cudaGetLastError cudaGetErrorString cudaPeekAtLastError cudaStreamCreate cudaStreamDestroy cudaStreamQuery cudaStreamSynchronize cudaStreamWaitEvent cudaEventCreate cudaEventCreateWithFlags cudaEventDestroy cudaEventElapsedTime cudaEventQuery cudaEventRecord cudaEventSynchronize
  syn keyword fortranIntrinsicCUDA    cudaFuncGetAttributes cudaFuncSetCacheConfig cudaFuncSetSharedMemConfig cudaSetDoubleForDevice cudaSetDoubleForHost cudaFree cudaFreeArray cudaFreeHost cudaGetSymbolAddress cudaGetSymbolSize
  syn keyword fortranIntrinsicCUDA    cudaHostAlloc cudaHostGetDevicePointer cudaHostGetFlags cudaHostRegister cudaHostUnregister cudaMalloc cudaMallocArray cudaMallocHost cudaMallocPitch cudaMalloc3D cudaMalloc3DArray
  syn keyword fortranIntrinsicCUDA    cudaMemcpy cudaMemcpyArraytoArray cudaMemcpyAsync cudaMemcpyFromArray cudaMemcpyFromSymbol cudaMemcpyFromSymbolAsync cudaMemcpyPeer cudaMemcpyPeerAsync cudaMemcpyToArray cudaMemcpyToSymbol cudaMemcpyToSymbolAsync cudaMemcpy2D cudaMemcpy2DArrayToArray cudaMemcpy2DAsync cudaMemcpy2DFromArray cudaMemcpy2DToArray cudaMemcpy3D cudaMemcpy3DAsync
  syn keyword fortranIntrinsicCUDA    cudaMemGetInfo cudaMemset cudaMemset2D cudaMemset3D cudaDeviceCanAccessPeer cudaDeviceDisablePeerAccess cudaDeviceEnablePeerAccess cudaPointerGetAttributes cudaDriverGetVersion cudaRuntimeGetVersion

  syn region none matchgroup=fortranType start="<<<" end=">>>" contains=ALLBUT,none
endif

syn cluster fortranCommentGroup contains=fortranTodo

syn match fortranContinueMark		display "&"

syn match fortranComment	excludenl "!.*$" contains=@fortranCommentGroup,@spell
syn match fortranOpenMP		excludenl 		"^\s*!\$\(OMP\)\=&\=\s.*$"

"cpp is often used with Fortran
syn match  cPreProc    "^\s*#\s*\(define\|ifdef\)\>.*"
syn match  cPreProc    "^\s*#\s*\(elif\|if\)\>.*"
syn match  cPreProc    "^\s*#\s*\(ifndef\|undef\)\>.*"
syn match  cPreProc    "__FILE__"
syn match  cPreProc    "__LINE__"
syn match  cPreCondit  "^\s*#\s*\(else\|endif\)\>.*"
syn region cIncluded   contained start=+"[^("]+ skip=+\\\\\|\\"+ end=+"+ contains=fortranLeftMargin,fortranContinueMark,fortranSerialNumber
"syn region	cIncluded	        contained start=+"[^("]+ skip=+\\\\\|\\"+ end=+"+
syn match	cIncluded		contained "<[^>]*>"
syn match	cInclude		"^\s*#\s*include\>\s*["<]" contains=cIncluded

"Synchronising limits assume that comment and continuation lines are not mixed
if exists("fortran_fold") || exists("fortran_more_precise")
  syn sync fromstart
else
  syn sync minlines=30
endif

if exists("fortran_fold")
  if exists("fortran_fold_multilinecomments")
    syn match fortranMultiLineComments transparent fold "\(^\s*!.*\(\n\|\%$\)\)\{4,}" contains=ALLBUT,fortranMultiCommentLines
  endif
endif

" Define the default highlighting.
" The default highlighting differs for each dialect.
" Transparent groups:
" fortranParen, fortranLeftMargin
" fortranProgram, fortranModule, fortranSubroutine, fortranFunction,
" fortranBlockData
" fortran77Loop, fortran90Loop, fortranIfBlock, fortranCase
" fortranMultiCommentLines
hi def link fortranKeyword         Keyword
hi def link fortranConstructName   Identifier
hi def link fortranConditional     Conditional
hi def link fortranRepeat          Repeat
hi def link fortranTodo            Todo
hi def link fortranContinueMark    Special
hi def link fortranString          String
hi def link fortranNumber          Number
hi def link fortranBinary          Number
hi def link fortranOctal           Number
hi def link fortranHex             Number
hi def link fortranOperator        Operator
hi def link fortranBoolean         Boolean
hi def link fortranLabelError      Error
hi def link fortranObsolete        Todo
hi def link fortranType            Type
hi def link fortranStructure       Type
hi def link fortranStorageClass    StorageClass
"=====================================================================[T-Flows]=
" Got rid of the special group for call ...
" hi def link fortranCall		Function
"---------------------------------------------------------------------[T-Flows]-
hi def link fortranUnitHeader      fortranPreCondit
hi def link fortranReadWrite       Keyword
hi def link fortranIO              Keyword
hi def link fortranIntrinsic       Function
hi def link fortranConstant        Constant

" To stop deleted & obsolescent features being highlighted as Todo items,
" comment out the next 5 lines and uncomment the 5 lines after that
hi def link fortranUnitHeaderOb    fortranObsolete
hi def link fortranKeywordOb       fortranObsolete
hi def link fortranConditionalOb   fortranObsolete
hi def link fortranTypeOb          fortranObsolete
hi def link fortranKeywordDel      fortranObsolete
"hi def link fortranUnitHeaderOb    fortranUnitHeader
"hi def link fortranKeywordOb       fortranKeyword
"hi def link fortranConditionalOb   fortranConditional
"hi def link fortranTypeOb          fortranType
"hi def link fortranKeywordDel      fortranKeyword

if b:fortran_dialect == "F"
  hi! def link fortranIntrinsicR         fortranObsolete
  hi! def link fortranUnitHeaderR        fortranObsolete
  hi! def link fortranTypeR              fortranObsolete
  hi! def link fortranStorageClassR      fortranObsolete
  hi! def link fortranOperatorR          fortranObsolete
  hi! def link fortranInclude            fortranObsolete
  hi! def link fortranLabelNumber        fortranObsolete
  hi! def link fortranTarget             fortranObsolete
  hi! def link fortranFloatIll           fortranObsolete
  hi! def link fortranIOR                fortranObsolete
  hi! def link fortranKeywordR           fortranObsolete
  hi! def link fortranStringR            fortranObsolete
  hi! def link fortranConditionalR       fortranObsolete
else
  hi! def link fortranIntrinsicR         fortranIntrinsic
  hi! def link fortranUnitHeaderR        fortranPreCondit
  hi! def link fortranTypeR              fortranType
  hi! def link fortranStorageClassR      fortranStorageClass
  hi! def link fortranOperatorR          fortranOperator
  hi! def link fortranInclude            Include
  hi! def link fortranLabelNumber        Special
  hi! def link fortranTarget             Special
  hi! def link fortranFloatIll           fortranFloat
  hi! def link fortranIOR                fortranIO
  hi! def link fortranKeywordR           fortranKeyword
  hi! def link fortranStringR            fortranString
  hi! def link fortranConditionalR       fortranConditional
endif

"=====================================================================[T-Flows]=
" T-Flows specific
hi def link fortranConstantTflows  fortranConstant
hi def link fortranTypeTflows      fortranType
hi def link fortranObjectTflows    fortranIntrinsic
hi def link fortranMpiTflows       fortranIntrinsic
hi def link fortranPetscTflows     fortranIntrinsic
hi def link fortranMacroTflows     fortranIntrinsic
hi def link fortranGlobalTflows    Todo
"---------------------------------------------------------------------[T-Flows]-

" CUDA
hi def link fortranIntrinsicCUDA   fortranIntrinsic
hi def link fortranTypeCUDA        fortranType
hi def link fortranStringCUDA      fortranString

hi def link fortranFormatSpec      Identifier
hi def link fortranFloat           Float
hi def link fortranPreCondit       PreCondit
hi def link cIncluded              fortranString
hi def link cInclude               Include
hi def link cPreProc               PreProc
hi def link cPreCondit             PreCondit
hi def link fortranOpenMP          PreProc
hi def link fortranParenError      Error
hi def link fortranComment         Comment
hi def link fortranSerialNumber    Todo
hi def link fortranTab             Error

" Uncomment the next line if you use extra intrinsics provided by vendors
"hi def link fortranExtraIntrinsic	Function

let b:current_syntax = "fortran"

let &cpo = s:cpo_save
unlet s:cpo_save
" vim: ts=8 tw=132
